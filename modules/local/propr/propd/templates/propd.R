#!/usr/bin/env Rscript


################################################
################################################
## Functions                                  ##
################################################
################################################

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse
parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Boolean. TRUE if first row is header. False without header.
#' @param row.names The first column is used as row names by default.
#' Otherwise, give another number. Or use NULL when no row.names are present.
#'
#' @return output Data frame
read_delim_flexible <- function(file, header = TRUE, row.names = 1, check.names = TRUE){

    ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1))

    if (ext == "tsv" || ext == "txt") {
        separator <- "\\t"
    } else if (ext == "csv") {
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }

    mat <- read.delim(
        file,
        sep = separator,
        header = header,
        row.names = row.names,
        check.names = check.names
    )
}

#' Get hub genes from adjacency matrix
#' 
#' Here hub genes are those that have a degree higher than the expected degree.
#' The expected degree is the number of connections that each gene would have
#' if the connections were distributed uniformly. In other words, the average 
#' degree by node.
#' 
#' @param adj Adjacency matrix
#' 
#' @return data frame with hub genes
get_hub_genes_from_adjacency <- function(adj){

    # get the expected degree
    degree_per_gene <- rowSums(adj)
    total_degree <- sum(degree_per_gene)
    n_nodes <- sum(degree_per_gene > 0)
    expected_degree <- total_degree / n_nodes

    # get hub genes
    hub_genes = degree_per_gene[degree_per_gene > expected_degree]
    hub_genes = data.frame(
        'feature' = names(hub_genes),
        'degree' = as.numeric(hub_genes)
    )
    names(hub_genes) <- c(opt\$features_id_col, 'degree')
    hub_genes <- hub_genes[order(hub_genes\$degree, decreasing=TRUE),]

    return(hub_genes)
}

################################################
################################################
## Parse arguments                            ##
################################################
################################################

# Set defaults and classes

opt <- list(
    prefix            = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    count             = '$count',
    samplesheet       = '$samplesheet',
    features_id_col   = 'gene_id',            # column name of feature ids
    obs_id_col        = 'sample',             # column name of observation ids
    group_col         = 'treatment',          # column name of grouping variable
    alpha             = NA,                   # alpha for boxcox transformation
    moderated         = TRUE,                 # use moderated theta
    fdr               = 0.05,                 # FDR threshold
    permutation       = 0,                    # if permutation > 0, use permutation test to compute FDR
    number_of_cutoffs = 100,                  # number of cutoffs for permutation test
    seed              = NA,                   # seed for reproducibility
    ncores            = as.integer('$task.cpus')
)
opt_types <- list(
    prefix            = 'character',
    count             = 'character',
    samplesheet       = 'character',
    features_id_col   = 'character',
    obs_id_col        = 'character',
    group_col         = 'character',
    alpha             = 'numeric',
    moderated         = 'logical',
    fdr               = 'numeric',
    permutation       = 'numeric',
    number_of_cutoffs = 'numeric',
    seed              = 'numeric',
    ncores            = 'numeric'
)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    } else {

        # Preserve classes from defaults where possible
        args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])

        # handle NA, and avoid errors when NA is provided by user as character
        if (args_opt[[ao]] %in% c('NA', NA)) args_opt[[ao]] <- NA

        # replace values
        opt[[ao]] <- args_opt[[ao]]
    }
}

# Check if required parameters have been provided

required_opts <- c('count','samplesheet')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]
if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid

for (file_input in c('count','samplesheet')){
    if (is.null(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }
    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

# TODO maybe add a function to pretty print the arguments?
print(opt)

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(propr)

################################################
################################################
## Perform differential proportionality       ##
################################################
################################################

# set seed when required

if (!is.na(opt\$seed)) {
    warning('Setting seed ', opt\$seed, ' for reproducibility')
    set.seed(opt\$seed)
}

# read matrix

mat <- read_delim_flexible(
    opt\$count,
    header = TRUE,
    row.names = opt\$features_id_col,
    check.names = FALSE
)
mat <- t(mat)  # transpose matrix to have features (genes) as columns

# parse group
# This creates a vector referring to the group id for each observation.
# The vector should have 2+ different groups, so that differential proportionality will 
# be computed to compare the variances between and within groups. TODO one can parse the 
# 'group_col' from the contrast file information as the other modules

samplesheet <- read_delim_flexible(
    opt\$samplesheet,
    header = TRUE,
    row.names = NULL,
    check.names = FALSE
)
tmp <- samplesheet[[opt\$group_col]]
names(tmp) <- samplesheet[[opt\$obs_id_col]]
group <- as.vector(tmp[rownames(mat)])
if (length(group) != nrow(mat)) stop('Error when parsing group')

# compute differential proportionality

pd <- propd(
    mat,
    group    = group,
    alpha    = opt\$alpha,
    weighted = FALSE,
    p        = opt\$permutation
)

# use F-stat FDR-adjusted p-values to get significant pairs, if permutation == 0
# otherwise, get FDR values using permutation tests (more computationally expensive but likely more conservative FDRs)

if (opt\$permutation == 0) {

    warning('FDR-adjusted p-values are used to get significant pairs.')

    # update FDR-adjusted p-values

    pd <- updateF(
        pd,
        moderated = opt\$moderated
    )
    if (opt\$moderated) pd <- setActive(pd, what='theta_mod')

    # get adjacency matrix

    adj <- getAdjacencyFstat(
        pd,
        pval=opt\$fdr,
        fdr_adjusted=TRUE
    )

    # get hub genes

    hub_genes <- get_hub_genes_from_adjacency(adj)

    # get significant pairs and classify them into red/yellow/green pairs

    results <- getSignificantResultsFstat(
        pd,
        pval=opt\$fdr,
        fdr_adjusted=TRUE
    )
    results <- results[,c("Partner", "Pair", "theta")]
    results\$class <- "red"
    results\$class[which(results\$Pair %in% hub_genes\$gene | results\$Partner %in% hub_genes\$gene)] <- "yellow"
    results\$class[which(results\$Pair %in% hub_genes\$gene & results\$Partner %in% hub_genes\$gene)] <- "green"

} else {

    warning('Permutation tests are used to compute FDR values.')

    # update FDR values using permutation tests

    pd <- updateCutoffs(
        pd,
        number_of_cutoffs = 100,
        ncores = opt\$ncores
    )

    # TODO take top n pairs when no cutoff has FDR below desired threshold
    cutoff <- getCutoffFDR(
        pd, 
        fdr=opt\$fdr,
        window_size=1
    )
    if (!cutoff) stop('No cutoff has FDR below desired threshold')

    # get adjacency matrix

    adj <- getAdjacencyFDR(
        pd,
        fdr=opt\$fdr,
        window_size=1
    )

    # get hub genes

    hub_genes <- get_hub_genes_from_adjacency(adj)

    # get significant pairs and classify them into red/yellow/green pairs

    results <- getSignificantResultsFDR(
        pd,
        fdr=opt\$fdr,
        window_size=1
    )
    results <- results[,c("Partner", "Pair", "theta")]
    results\$class <- "red"
    results\$class[which(results\$Pair %in% hub_genes\$gene | results\$Partner %in% hub_genes\$gene)] <- "yellow"
    results\$class[which(results\$Pair %in% hub_genes\$gene & results\$Partner %in% hub_genes\$gene)] <- "green"
}

################################################
################################################
## Generate outputs                           ##
################################################
################################################

saveRDS(
    pd,
    file = paste0(opt\$prefix, '.propd.rds')
)

write.table(
    getResults(pd),
    file      = paste0(opt\$prefix, '.propd.results.tsv'),
    col.names = TRUE,
    row.names = FALSE,
    sep       = '\\t',
    quote     = FALSE
)

write.table(
    results,
    file      = paste0(opt\$prefix, '.propd.results_filtered.tsv'),
    col.names = TRUE,
    row.names = FALSE,
    sep       = '\\t',
    quote     = FALSE
)

write.table(
    adj,
    file      = paste0(opt\$prefix, '.propd.adjacency.csv'),
    col.names = TRUE,
    row.names = TRUE,
    sep       = ',',
    quote     = FALSE
)

write.table(
    hub_genes,
    file      = paste0(opt\$prefix, '.propd.hub_genes.tsv'),
    col.names = TRUE,
    row.names = FALSE,
    sep       = '\\t',
    quote     = FALSE
)

if (opt\$permutation > 0) {
    write.table(
        pd@fdr,
        file      = paste0(opt\$prefix, '.propd.fdr.tsv'),
        col.names = TRUE,
        sep       = '\\t',
        quote     = FALSE
    )
}

################################################
################################################
## WARNINGS                                   ##
################################################
################################################

sink(paste0(opt\$prefix, ".warnings.log"))
print(warnings())
sink()

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste0(opt\$prefix, ".R_sessionInfo.log"))
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

propr.version <- as.character(packageVersion('propr'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-propr:', propr.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################