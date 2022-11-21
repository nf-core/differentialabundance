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
    args_vals <- unlist(lapply(args_list, function(y) strsplit(y, ' +')))

    as.list(structure(args_vals[c(FALSE, TRUE)], names = args_vals[c(TRUE, FALSE)]))
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
    gene_id_col = 'gene_id',
    sample_file = '$samplesheet',
    sample_id_col = 'experiment_accession'
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
sample.sheet <- read.table(file = opt\$sample_file, sep="\t", header=T)
count.table <- read.table(file = opt\$count_file, sep="\t", header=T)

# Rewrite input tables to a standard differential abundance input; columns in sample sheet need to be renamed from file names to sample names
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

# Write modified count table
write.table(count.table,
    file = 'modified_counts.tsv',
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE)

# Write samplesheet as CSV
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
# TODO
output_prefix <- "dummy_prefix"
sink(paste(output_prefix, "R_sessionInfo.log", sep = '.'))
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
