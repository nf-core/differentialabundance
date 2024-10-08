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

    return(mat)
}

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
    features_id_col   = 'gene_id',                  # column name for features (genes)

    # parameters for computing correlation coefficients
    metric            = 'rho',                      # correlation metric: rho, phi, phs, cor, vlr, pcor, pcor.shrink, pcor.bshrink
    ivar              = 'clr',                      # transformation: clr, alr, or the name(s) or index(es) of the variable(s) to be used as reference(s)
    alpha             = NA,                         # alpha value for Box-Cox transformation

    # parameters for getting the significant coefficients
    fdr               = 0.05,                       # FDR threshold
    tails             = "right",                    # FDR tail: right or both
    permutation       = 100,                        # number of permutations
    number_of_cutoffs = 100,                        # number of cutoffs for which compute the FDRs

    # other parameters
    seed              = NA,                         # seed for reproducibility
    ncores            = as.integer('$task.cpus')
)

opt_types <- list(
    prefix            = 'character',
    count             = 'character',
    features_id_col   = 'character',
    metric            = 'character',
    ivar              = 'character',
    alpha             = 'numeric',
    fdr               = 'numeric',
    tails             = 'character',
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

        # Preserve classes from defaults
        args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])

        # handle NA, and avoid errors when NA is provided by user as character
        if (args_opt[[ao]] %in% c('NA', NA)) args_opt[[ao]] <- NA

        # replace values
        opt[[ao]] <- args_opt[[ao]]
    }
}

# Check if required parameters have been provided

required_opts <- c('count')  # only count data is strictly required, other parameters have defaults
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]
if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid

for (file_input in c('count')){
    if (is.null(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }
    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

# check parameters

if (!opt\$metric %in% c('rho', 'phi', 'phs', 'cor', 'vlr', 'pcor', 'pcor.shrink', 'pcor.bshrink')) {
    stop('Please make sure you provided the correct metric')
}

if (opt\$metric == 'pcor.bshrink'){
    if (!is.na(opt\$alpha)) stop('Box-cox transformation is not implemented for pcor.bshrink yet.')
    if (!opt\$ivar %in% c('clr', 'alr')) stop('Please make sure you provided the correct transformation: clr or alr')

} else {
    if (is.na(opt\$ivar)) warning('Warning: No transformation is required by user. We assume the input count data was already properly transformed.')
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
## Perform correlation analysis               ##
################################################
################################################

# set seed when required

if (!is.na(opt\$seed)) {
    warning('Setting seed ', opt\$seed, ' for reproducibility')
    set.seed(opt\$seed)
}

# load count matrix

mat <- read_delim_flexible(
    opt\$count,
    header = TRUE,
    row.names = opt\$features_id_col,
    check.names = FALSE
)
mat <- t(mat)  # transpose matrix to have features (genes) as columns

# Compute correlation coefficients

pr <- propr(
    mat,
    metric  = opt\$metric,
    ivar    = opt\$ivar,
    alpha   = opt\$alpha,
    p       = opt\$permutation
)

if (opt\$permutation > 0) {

    # update FDRs for each coefficient cutoff

    pr <- updateCutoffs(
        pr,
        number_of_cutoffs=opt\$number_of_cutoffs,
        tails=opt\$tails,
        ncores=opt\$ncores
    )

    # get cutoff at given FDR threshold

    cutoff <- getCutoffFDR(
        pr,
        fdr=opt\$fdr,
        window_size=1
    )

    if (cutoff) {

        # get adjacency matrix with the significant edges

        adj <- getAdjacencyFDR(
            pr,
            fdr=opt\$fdr,
            window_size=1
        )

    } else {
        # TODO take top n pairs when no cutoff has FDR below desired threshold
        # For the moment, we just print a warning and set adj to NULL
        warning('No significant results found at FDR threshold ', opt\$fdr)
        adj <- NULL
    }
}

################################################
################################################
## Generate outputs                           ##
################################################
################################################

saveRDS(
    pr,
    file = paste0(opt\$prefix, '.propr.rds')
)

write.table(
    round(pr@matrix, 8),  # round matrix decimals to avoid floating point inconsistencies
    file      = paste0(opt\$prefix, '.propr.matrix.csv'),
    col.names = TRUE,
    row.names = TRUE,
    sep       = ',',
    quote     = FALSE
)

if (!is.null(adj)) {
    write.table(
        adj,
        file      = paste0(opt\$prefix, '.propr.adjacency.csv'),
        col.names = TRUE,
        row.names = FALSE,
        sep       = '\t',
        quote     = FALSE
    )
}

if (opt\$permutation > 0) {
    write.table(
        pr@fdr,
        file      = paste0(opt\$prefix, '.propr.fdr.tsv'),
        col.names = TRUE,
        row.names = FALSE,
        sep       = '\t',
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
