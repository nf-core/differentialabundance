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
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#'
#' @return output Data frame
read_delim_flexible <- function(file, header = TRUE, row.names = NULL){
    ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1))
    if (ext == "tsv" || ext == "txt") {
        separator <- "\\t"
    } else if (ext == "csv") {
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }
    read.delim(
        file,
        sep = separator,
        header = header,
        row.names = row.names
    )
}
#' Round numeric dataframe columns to fixed decimal places by applying
#' formatting and converting back to numerics
#'
#' @param dataframe A data frame
#' @param columns Which columns to round (assumes all of them by default)
#' @param digits How many decimal places to round to?
#'
#' @return output Data frame
round_dataframe_columns <- function(df, columns = NULL, digits = 8){
    if (is.null(columns)){
        columns <- colnames(df)
    }
    df[,columns] <- format(data.frame(df[, columns]), nsmall = digits)
    # Convert columns back to numeric
    for (c in columns) {
        df[[c]][grep("^ *NA\$", df[[c]])] <- NA
        df[[c]] <- as.numeric(df[[c]])
    }
    df
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

# I've defined these in a single array like this so that we could go back to an
# optparse-driven method in future with module bin/ directories, rather than
# the template

# Set defaults and classes

opt <- list(
    count_file = '$counts',
    gene_id_col = 'Geneid',
    sample_file = '$samplesheet',
    sample_id_col = 'Sample Name'
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }else{
        opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
    }
}

# Check file inputs are valid

for (file_input in c('count_file', 'sample_file')){
    if (is.null(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }

    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

################################################
################################################
## READ IN COUNTS FILE AND SAMPLE METADATA    ##
################################################
################################################
count.table <- read.table(file = opt\$count_file, sep="\t", header=T, check.names=F)
sample.sheet <- read.table(file = opt\$sample_file, sep="\t", header=T, check.names=F)

# Check if column names are found in tables

if (!(opt\$gene_id_col %in% names(count.table))) {
    stop(paste("Could not find gene_id_col", opt\$gene_id_col, "in the following count matrix columns:\n",
                toString(names(count.table))))
}
if (!(opt\$sample_id_col %in% names(sample.sheet))) {
    stop(paste("Could not find sample_id_col", opt\$sample_id_col, "in the following input sample sheet columns:\n",
                toString(names(sample.sheet))))
}

# Rewrite input tables to a standard differential abundance input; columns in sample sheet need to be
# renamed from file names to sample names

col_names <- list(colnames(count.table))[[1]]
remaining_names <- col_names
sample_names <- sample.sheet[opt\$sample_id_col][[1]]

# Order sample names by decreasing length so that longer names are checked first;
# this prevents Sample11.Aligned, Sample.12.Aligned etc. all being renamed to Sample1

sample_names <- sample_names[order(nchar(sample_names), decreasing=T)]
for (c in c(1:length(col_names))) {
    for (q in c(1:length(sample_names))) {
        if (grepl(toString(sample_names[q]), col_names[[c]], fixed = TRUE)) {
            old_name <- col_names[[c]]
            new_name <- toString(sample_names[q])
            names(count.table)[names(count.table) == old_name] <- new_name
            remaining_names <- remaining_names[remaining_names != old_name]
            break
        }
    }
}

# Remove non-sample columns, then check if all samples were renamed

remaining_names <- remaining_names[remaining_names != opt\$gene_id_col]
remaining_names <- remaining_names[remaining_names != "gene_name"]
if (length(remaining_names > 0)) {
    print(paste0("Warning: Could not rename the following columns:\n", toString(remaining_names)))
}

# Write modified count table as TSV

write.table(count.table,
    file = 'modified_counts.tsv',
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE)

# Write modified samplesheet as CSV

write.table(sample.sheet,
    file = 'modified_samplesheet.csv',
    col.names = TRUE,
    row.names = FALSE,
    sep = ',',
    quote = FALSE)

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink("R_sessionInfo.log")
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
