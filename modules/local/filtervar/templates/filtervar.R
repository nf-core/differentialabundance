

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

    if ( (row.names == 'gene_id') & ('gene_name' %in% colnames(mat)) ){
        mat <- mat[, -which(colnames(mat) == 'gene_name')]
    } else if ( (row.names == 'gene_name') & ('gene_id' %in% colnames(mat)) ){
        mat <- mat[, -which(colnames(mat) == 'gene_id')]
    }

    return(mat)
}

read_delim_flexible2 <- function(file, header = TRUE){

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
        header = header
    )
    return(mat)
}



################################################
################################################
## Parse arguments                            ##
################################################
################################################

opt <- list(
    count            = '$count',
    prefix           = ifelse('$task.ext.prefix' == 'null', '$meta.pathway_name', '$task.ext.prefix'),
    transformation   = 'clr',
    reference        = NA,
    alpha            = NA,
    metric           = 'pcor.bshrink',
    permutation      = 0,
    cutoff_min       = NA,
    cutoff_max       = NA,
    cutoff_interval  = NA,
    ncores           = as.integer('$task.cpus'),
    features_id_col  = 'gene_id',
    fixseed          = FALSE,
    adjacency        = FALSE,
    fdrVal           = 0.05,
    adj_matrix       = '$adj_matrix',
    filterVar        = 'yes'
)
opt_types <- list(
    count            = 'character',
    prefix           = 'character',
    transformation   = 'character',
    reference        = 'character',
    alpha            = 'numeric',
    metric           = 'character',
    permutation      = 'numeric',
    cutoff_min       = 'numeric',
    cutoff_max       = 'numeric',
    cutoff_interval  = 'numeric',
    ncores           = 'numeric',
    features_id_col  = 'character',
    fixseed          = 'logical',
    adjacency        = 'logical',
    fdrVal           = 'numeric',
    adj_matrix       = 'character',
    filterVar        = 'character'
)


# Apply parameter overrides
args_opt <- parse_args('$task.ext.args')

for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    } else {

        # Preserve classes from defaults where possible
        if (! is.null(opt[[ao]])){
            args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        }
        # set NA
        if (args_opt[[ao]] %in% c('NA', NA, 'null')){
            args_opt[[ao]] <- NA
        }
        opt[[ao]] <- args_opt[[ao]]
    }
}

# Check if required parameters have been provided
required_opts <- c('count')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]
if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

################################################
################################################
## Perform variable selection                 ##
################################################
################################################

# read matrix
A <- read_delim_flexible(
    opt\$adj_matrix,
    header = TRUE,
    row.names = 1,
    check.names = TRUE
)

count <- read_delim_flexible2(
    opt\$count,
    header = TRUE
)

### Determine most differentially proportional genes

# Set diagonal in A to 0
diag(A) <- 0

# Sum values in adjacency and add as an extra column
per_gene_connection <- rowSums(A)

A\$per_gene <- per_gene_connection

A <- A[order(A\$per_gene, decreasing = TRUE),]

# Define selection criteria 

max_gene_number <- ncol(count)*10 # 10x samples for technical reasons (pcor)

#Calculate connection threshold 
total_connections <- sum(per_gene_connection)/2 # 2 because the matrix is symmetric
possible_connections <- nrow(count)*(nrow(count)-1)/2

percentage_expected <- total_connections/possible_connections
connection_threshold <- percentage_expected * nrow(count)

# Filter count matrix according to selected genes 

col_genes <- which(names(count) == opt\$features_id_col)

if (opt\$filterVar == 'yes'){
    # select only differentially proportional genes 
    top_genes <- rownames(A[which(A\$per_gene > connection_threshold),])
    count_filtered <- count[count[,col_genes] %in% top_genes,]
    warning("non differentially proportional genes were removed before correlation analysis")

} else if (max_gene_number < nrow(count) & opt\$metric== 'pcor.bshrink'){
    # select the maximum number of genes to perform partial correlation
    top_genes <- rownames(A[1:gene_number,])
    count_filtered <- count[count[,col_genes] %in% top_genes,]
    warning("some genes were removed to perform partial correlation")

}else{
    # no genes were removed
    count_filtered <- count
    warning("No genes were removed")
}


################################################
################################################
## Generate outputs                           ##
################################################
################################################

write.table(
    count_filtered,
    file      = paste0(opt\$prefix, '.count_filtered.tsv'),
    col.names = TRUE,
    row.names = FALSE,
    sep       = '\t',
    quote     = FALSE
)

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

