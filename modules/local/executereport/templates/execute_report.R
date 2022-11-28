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

opt <- list(
    report_file = '$report_file',
    args = '$args'
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

for (file_input in c('report_file')){
    if (is.null(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }

    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(rmarkdown)

###########

wd=getwd()
rmarkdown::render(opt\$report_file, output_file = "differentialabundance_report.html", params = list(args = opt\$args), knit_root_dir = wd, output_dir = wd)
file.copy("differentialabundance_report.html", "/home-link/iivow01/git/differentialabundance/error/differentialabundance_report.html")
