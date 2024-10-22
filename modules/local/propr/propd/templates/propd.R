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

#' Get connectivity of genes from adjacency matrix
#'
#' The connectivity of a gene is the number of connections it has with other genes.
#' In other words, the degree of a gene. This connectivity can be weighted. In this
#' way, the strength of the theta value can be taken into account.
#' These information can be used to summarize the network and identify genes that are
#' potentially changing between groups.
#'
#' @param pd propd object
#' @param adj Adjacency matrix
#' @param de Differential proportionality values for each gene with respect to the
#' normalization reference
#' @param cutoff Cutoff value for de values to be considered significant
#'
#' @return data frame with the following columns: gene_id, degree, weighted_degree,
#' genewise_theta, average_theta, classification
get_connectivity <- function(pd, adj, de, cutoff, features_id_col='gene_id'){

    # initialize empty data frame
    connectivity <- data.frame(matrix(NA, nrow=ncol(pd@counts), ncol=6))
    colnames(connectivity) <- c(
        features_id_col,
        'degree',
        'weighted_degree',
        'genewise_theta',
        'average_theta',
        'classification'
    )

    # add features ids
    connectivity[,1] <- colnames(pd@counts)

    # add degree
    # degree is the number of connections a gene has with other genes
    connectivity[,2] <- colSums(adj)

    # add weighted degree
    # each connection is weighted by the theta value
    # so lower theta values (higher between group variance than within group variance) will have a higher weight
    # NOTE this is a placeholder for the proper weighted degree, maybe we are gonna change the way how we compute it
    mat <- getMatrix(pd)
    diag(mat) <- NA
    connectivity[,3] <- colSums((1 - mat) * adj, na.rm=TRUE)

    # add genewise theta
    # a theta value can be associated to each gene by calculating the between group variance vs within group variance
    # of the gene normalized with respect to a reference (in this case the geometric mean of the sample)
    connectivity[,4] <- de

    # add average theta of the connections
    connectivity[,5] <- colSums(mat * adj, na.rm=TRUE) / colSums(adj)

    # classification
    # green for DE genes, and red for non-DE genes
    connectivity[,6] <- 'green'
    connectivity[which(de > cutoff), 6] <- 'red'

    return(connectivity)
}

#' Determine hub genes based on connectivity
#'
#' Here hub genes are those that have a degree higher than the expected degree.
#' The expected degree is the number of connections that each gene would have
#' if the connections were distributed uniformly. In other words, the average
#' degree by node.
#'
#' @param connectivity Data frame with connectivity
<<<<<<< Updated upstream
#' @param cutoff Theta value for which DP pairs are considered significant.
#' @param weighted Boolean. If TRUE, use weighted degree to determine hub genes.
#' Otherwise, use degree.
#'
#' @return filtered and sorted connectivity data frame with hub genes
get_hub_genes <- function(connectivity, cutoff, weighted=FALSE){

    # get the expected degree
    total_degree <- if (weighted) sum(connectivity\$weighted_degree) else sum(connectivity\$degree)
    n_nodes <- sum(connectivity\$degree > 0)
=======
#'
#' @return filtered connectivity data frame with hub genes
get_hub_genes <- function(connectivity){

    # get the expected degree
    total_degree <- sum(connectivity\$degree)
    n_nodes <- sum(connectivity > 0)
>>>>>>> Stashed changes
    expected_degree <- total_degree / n_nodes

    # get hub genes
    hub_genes <- connectivity[which(connectivity\$degree > expected_degree),]

<<<<<<< Updated upstream
    # sort hub genes
    if (weighted) {
        hub_genes <- hub_genes[order(hub_genes\$weighted_degree, decreasing=TRUE),]
    } else {
        hub_genes <- hub_genes[order(hub_genes\$degree, decreasing=TRUE),]
    }

    return(hub_genes)
}

#' Plot pairs of genes
#'
#' This function plots the following pairs of genes:
#' - xy vs sample
#' - x vs y
#' - x vs sample
#' - y vs sample
#' - xr vs yr
#' - xr vs sample
#' - yr vs sample
#' The pairs are colored according to the group they belong to.
#'
#' @param df Data frame with the following columns: xy, x, y, xr, yr, sample, group, color
#' @param x Name of the gene x
#' @param y Name of the gene y
#' @param title Title of the plot
plot_pairs <- function(df, x, y, title){

    # Define the layout
    layout_matrix <- matrix(c(
        0, 1, 8,    # First row
        2, 3, 4,    # Second row
        5, 6, 7     # Third row
        ), nrow = 3, ncol = 3, byrow = TRUE)
    layout(layout_matrix, widths = c(1, 1, 1), heights = c(1, 1, 1))

    # Adjust margins and text sizes
    par(mar = c(4, 4, 1, 1), cex = 1.5, lwd = 1.5)

    # plot xy vs sample
    plot(x=df\$sample, y=df\$xy, xlab='sample', ylab=paste0(x, '/', y), col=df\$color)

    # plot x vs y
    plot(x=df\$x, y=df\$y, xlab=x, ylab=y, col=df\$color)

    # plot x vs sample
    plot(x=df\$sample, y=df\$x, xlab='sample', ylab=x, col=df\$color)

    # plot y vs sample
    plot(x=df\$sample, y=df\$y, xlab='sample', ylab=y, col=df\$color)

    # plot xr vs yr
    plot(x=df\$xr, y=df\$yr, xlab=paste0(x, '/ref'), ylab=paste0(y, '/ref'), col=df\$color)

    # plot xr vs sample
    plot(x=df\$sample, y=df\$xr, xlab='sample', ylab=paste0(x, '/ref'), col=df\$color)

    # plot yr vs sample
    plot(x=df\$sample, y=df\$yr, xlab='sample', ylab=paste0(y, '/ref'), col=df\$color)

    # add legend
    plot.new()
    par(mar = c(0, 0, 0, 0))
    legend(
        "center",
        legend = unique(df\$group),
        col = unique(df\$color),
        pch = 19,
        cex = 1.5,
        bty = "n")

    # TODO the title does not appear somehow, fix it
    # Add main title
    mtext(
        title,
        side = 3,
        outer = TRUE,
        line = 1,
        cex = 2,
        font = 2)
}

=======
    return(hub_genes)
}

>>>>>>> Stashed changes
################################################
################################################
## Parse arguments                            ##
################################################
################################################

# Set defaults and classes

opt <- list(
    prefix            = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),

    # input count matrix
    count             = '$count',
    features_id_col   = 'gene_id',            # column name of feature ids

    # comparison groups
    samplesheet       = '$samplesheet',
    obs_id_col        = 'sample',             # column name of observation ids
<<<<<<< Updated upstream
    contrast_variable = "$contrast_variable", # column name of contrast variable
    reference_group   = "$reference",         # reference group for contrast variable
    target_group      = "$target",            # target group for contrast variable
=======
    group_col         = 'treatment',          # column name of grouping variable
>>>>>>> Stashed changes

    # parameters for computing differential proportionality
    alpha             = NA,                   # alpha for boxcox transformation
    moderated         = TRUE,                 # use moderated theta

    # parameters for getting the significant differentially proportional pairs
    fdr               = 0.05,                 # FDR threshold
    permutation       = 0,                    # if permutation > 0, use permutation test to compute FDR
    number_of_cutoffs = 100,                  # number of cutoffs for permutation test

<<<<<<< Updated upstream
    # parameters for getting the hub genes
    weighted_degree   = FALSE,                 # use weighted degree for hub genes or not

=======
>>>>>>> Stashed changes
    # other parameters
    seed              = NA,                   # seed for reproducibility
    ncores            = as.integer('$task.cpus')
)

opt_types <- list(
    prefix            = 'character',
    count             = 'character',
    samplesheet       = 'character',
    features_id_col   = 'character',
    obs_id_col        = 'character',
<<<<<<< Updated upstream
    contrast_variable = 'character',
    reference_group   = 'character',
    target_group      = 'character',
=======
    group_col         = 'character',
>>>>>>> Stashed changes
    alpha             = 'numeric',
    moderated         = 'logical',
    fdr               = 'numeric',
    permutation       = 'numeric',
    number_of_cutoffs = 'numeric',
<<<<<<< Updated upstream
    weighted_degree   = 'logical',
=======
>>>>>>> Stashed changes
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

<<<<<<< Updated upstream
required_opts <- c('count','samplesheet','contrast_variable','reference_group','target_group')
=======
required_opts <- c('count','samplesheet')
>>>>>>> Stashed changes
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

<<<<<<< Updated upstream
# check parameters are valid

if (opt\$permutation < 0) {
    stop('permutation should be a positive integer')
}

=======
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
# and filter matrix and group values, so that only the contrasted groups are kept
# TODO propd can also handle more than two groups
# but that dont work properly with the contrast format
# Should we provide an alternative way to do that?
=======
# This creates a vector referring to the group id for each observation.
# The vector should have 2+ different groups, so that differential proportionality will
# be computed to compare the variances between and within groups. TODO one can parse the
# 'group_col' from the contrast file information as the other modules
>>>>>>> Stashed changes

samplesheet <- read_delim_flexible(
    opt\$samplesheet,
    header = TRUE,
    row.names = NULL,
    check.names = FALSE
)
<<<<<<< Updated upstream
samplesheet <- samplesheet[,c(opt\$obs_id_col, opt\$contrast_variable)]
idx <- which(samplesheet[,2] %in% c(opt\$reference_group, opt\$target_group))
mat <- mat[idx,]
samplesheet <- samplesheet[idx,]
group <- as.vector(samplesheet[,2])
if (length(group) != nrow(mat)) stop('Error when parsing group')
if (length(unique(group)) != 2) stop('Only two groups are allowed for contrast')
=======
tmp <- samplesheet[[opt\$group_col]]
names(tmp) <- samplesheet[[opt\$obs_id_col]]
group <- as.vector(tmp[rownames(mat)])
if (length(group) != nrow(mat)) stop('Error when parsing group')
>>>>>>> Stashed changes

# compute differential proportionality

pd <- propd(
    mat,
    group    = group,
    alpha    = opt\$alpha,
    weighted = FALSE,
    p        = opt\$permutation
)

<<<<<<< Updated upstream
# compute DE theta values
# this is the theta value for each gene with respect to the normalization reference
# in this case, the reference is the geometric mean of the sample
# These DE values are only for interpreting purposes.
# TODO if we want to use the outcome from other DE analysis, at some point we should
# divide the below part into maybe a separate module that can take the DE and DP values as input
# and coordinate them through the pipeline

ref <- exp(rowMeans(log(pd@counts)))
de <- runNormalization(pd, ref)

=======
>>>>>>> Stashed changes
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

<<<<<<< Updated upstream
    # get theta value for which FDR is below desired threshold
    # theta_cutoff is FALSE when no theta value has FDR below desired threshold
    # otherwise it is the theta value for which FDR is below desired threshold
    # Only when there is a meaningful theta, we can compute the next steps
    # that involve extracting the significant pairs.

    theta_cutoff <- getCutoffFstat(
=======
    # get adjacency matrix

    adj <- getAdjacencyFstat(
>>>>>>> Stashed changes
        pd,
        pval=opt\$fdr,
        fdr_adjusted=TRUE
    )
<<<<<<< Updated upstream
    if (theta_cutoff) {

        warning('Significant theta value found: ', theta_cutoff)

        # get adjacency matrix
        # this matrix will have 1s for significant pairs and 0s for the rest
        # diagonals are set to 0

        adj <- getAdjacencyFstat(
            pd,
            pval=opt\$fdr,
            fdr_adjusted=TRUE
        )

        # calculate gene connectivity and get hub genes

        connectivity <- get_connectivity(
            pd,
            adj,
            de,
            theta_cutoff,
            features_id_col=opt\$features_id_col
        )
        hub_genes <- get_hub_genes(connectivity, weighted=opt\$weighted_degree)

        # get significant pairs and classify them into red/yellow/green pairs

        results <- getSignificantResultsFstat(
            pd,
            pval=opt\$fdr,
            fdr_adjusted=TRUE
        )
        results <- results[,c("Partner", "Pair", "theta")]
        results\$classification <- "red"
        results\$classification[which(de[results\$Pair] < theta_cutoff | de[results\$Partner] < theta_cutoff)] <- "yellow"
        results\$classification[which(de[results\$Pair] < theta_cutoff & de[results\$Partner] < theta_cutoff)] <- "green"

        # sort significant pairs

        results <- results[order(results\$theta),]
    }
=======

    # calculate gene connectivity and get hub genes

    connectivity <- get_connectivity(adj)
    hub_genes <- get_hub_genes(connectivity)

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
>>>>>>> Stashed changes

} else {

    warning('Permutation tests are used to compute FDR values.')

    # update FDR values using permutation tests

    pd <- updateCutoffs(
        pd,
        number_of_cutoffs = 100,
        ncores = opt\$ncores
    )

<<<<<<< Updated upstream
    # get theta cutoff

    theta_cutoff <- getCutoffFDR(
=======
    # get cutoff
    # this is the cutoff used to get the significant pairs and ensemble of adjacency matrix

    cutoff <- getCutoffFDR(
>>>>>>> Stashed changes
        pd,
        fdr=opt\$fdr,
        window_size=1
    )
<<<<<<< Updated upstream
    if (theta_cutoff) {

        warning('Significant theta value found: ', theta_cutoff)
=======

    if (cutoff) {
>>>>>>> Stashed changes

        # get adjacency matrix

        adj <- getAdjacencyFDR(
            pd,
            fdr=opt\$fdr,
            window_size=1
        )

        # calculate gene connectivity and get hub genes

<<<<<<< Updated upstream
        connectivity <- get_connectivity(
            pd,
            adj,
            de,
            theta_cutoff,
            features_id_col=opt\$features_id_col
        )
        hub_genes <- get_hub_genes(connectivity, weighted=opt\$weighted_degree)
=======
        connectivity <- get_connectivity(adj)
        hub_genes <- get_hub_genes(connectivity)
>>>>>>> Stashed changes

        # get significant pairs and classify them into red/yellow/green pairs

        results <- getSignificantResultsFDR(
            pd,
            fdr=opt\$fdr,
            window_size=1
        )
        results <- results[,c("Partner", "Pair", "theta")]
<<<<<<< Updated upstream
        results\$classification <- "red"
        results\$classification[which(de[results\$Pair] < theta_cutoff | de[results\$Partner] < theta_cutoff)] <- "yellow"
        results\$classification[which(de[results\$Pair] < theta_cutoff & de[results\$Partner] < theta_cutoff)] <- "green"

        # sort significant pairs

        results <- results[order(results\$theta),]
    }
}

# deal with the situation when no significant thetas are found
# For the moment, we just print a warning and set adj, hub_genes and results to NULL
# TODO take top n pairs when no cutoff has FDR below desired threshold

if (!theta_cutoff) {
    warning('No theta value has FDR below desired threshold.')
    adj <- NULL
    connectivity <- NULL
    hub_genes <- NULL
    results <- NULL
}

=======
        results\$class <- "red"
        results\$class[which(results\$Pair %in% hub_genes\$gene | results\$Partner %in% hub_genes\$gene)] <- "yellow"
        results\$class[which(results\$Pair %in% hub_genes\$gene & results\$Partner %in% hub_genes\$gene)] <- "green"

    } else {
        # TODO take top n pairs when no cutoff has FDR below desired threshold
        # For the moment, we just print a warning and set adj, hub_genes and results to NULL
        warning('No pairs have FDR below desired threshold.')
        adj <- NULL
        connectivity <- NULL
        hub_genes <- NULL
        results <- NULL
    }
}

>>>>>>> Stashed changes
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

<<<<<<< Updated upstream
if (theta_cutoff) {
=======
if (!is.null(adj)) {
>>>>>>> Stashed changes
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
        connectivity,
        file      = paste0(opt\$prefix, '.propd.connectivity.tsv'),
        col.names = TRUE,
        row.names = FALSE,
        sep       = '\\t',
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
}

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
<<<<<<< Updated upstream
## Plot red pairs                             ##
################################################
################################################

if (theta_cutoff){

    # get ratios between each gene and the normalization reference
    ratios <- exp(logratio(pd@counts, 'clr', NA))

    # plot for each pair type
    for (pair_type in c("red", "yellow", "green")){
        pdf(paste0(opt\$prefix, '.propd.', pair_type, '_pairs.pdf'), width = 18, height = 18)

        # get pairs
        pairs <- results[results\$classification == pair_type,]

        # for the top pairs
        for (idx in c(1:3)){

            # get x and y genes
            x <- pairs[idx,'Partner']
            y <- pairs[idx,'Pair']

            # create data frame
            df <- data.frame(
                xy=pd@counts[,x]/pd@counts[,y],
                x=pd@counts[,x],
                y=pd@counts[,y],
                xr=ratios[,x],
                yr=ratios[,y],
                sample=c(1:nrow(pd@counts)),
                group=group,
                color=ifelse(group == opt\$target_group, 'red', 'blue'))
            df <- df[order(df\$group, df\$sample),]

            # plot
            title <- paste0("top ", idx, " ", pair_type, " pair with theta=", round(results[idx, 'theta'], 6))
            plot_pairs(df, x, y, title)
        }

    }
}

################################################
################################################
=======
>>>>>>> Stashed changes
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
