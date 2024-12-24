#!/usr/bin/env Rscript

# Call shinyngs parsing functions to validate simple matrix inputs

library(optparse)
library(yaml)

#Functions from accesory

#' Read tables of differential statistics
#'
#' @param filename File name of file with table of differential statistics
#' @param feature_id_column Column of stats file with feature identifiers
#' @param pval_column Column of stats file with p values
#' @param qval_column Column of stats file with adjust p values/ q values
#' @param fc_column Column of stats with fold changes
#' @param unlog_foldchanges Reverse a log on fold changes? Set to TRUE if values
#'   are logged.
#'
#' @return output Validated selected columns of differential stats files as a
#'   data frame
#' @export

read_differential_yaml <- function(filename,
                                                                                                                        feature_id_column = NULL,
                                                                                                                        pval_column = NULL,
                                                                                                                        qval_column = NULL,
                                                                                                                        fc_column = NULL,
                                                                                                                        unlog_foldchanges = FALSE) {
        st <- read_metadata(filename, id_col = feature_id_column)
        stats_cols <- c(feature_id_column, pval_column, qval_column, fc_column)

        success <-
                checkListIsSubset(
                        stats_cols,
                        colnames(st),
                        "stats variables",
                        "available stats columns"
                )
        st <- st[, stats_cols]

        for (c in colnames(st)) {
                st[[c]][grep("^ *NA$", st[[c]])] <- NA
                if (c %in% c(pval_column, qval_column, fc_column)) {
                        st[[c]] <- as.numeric(st[[c]])
                }
        }

        if (unlog_foldchanges) {
                st[[fc_column]] <- sign(st[[fc_column]]) * 2^(abs(st[[fc_column]]))
        }
        st
}

#' Take two delimiter-separated strings and generate a named vector
#'
#' @param elements_string String to be converted to vector elements
#' @param names_string String to be converted to vector names by default taken
#'   from elements_string.
#' @param sep Separator to use
#' @param prettify_names Boolean. Prettify element names?
#' @param simplify_files If elements are file and to be used to generate names,
#'   should we simplify by striping path and extension?
#'
#' @return output Named character vector
#' @export

stringsToNamedVector_yaml <- function(elements_string, names_string = NULL, sep = ",", prettify_names = TRUE, simplify_files = FALSE) {
        elements <- simpleSplit(elements_string, sep = sep)

        if (is.null(names_string)) {
                element_names <- elements
                if (simplify_files) {
                        element_names <- tools::file_path_sans_ext(gsub("_", " ", basename(element_names)))
                }
        } else {
                element_names <- simpleSplit(names_string, sep = sep)
        }

        if (length(elements) != length(element_names)) {
                stop(paste("List in", elements_string, "a different length to", names_string))
        }
        if (prettify_names) {
                element_names <- prettifyVariablename(element_names)
        }

        names(elements) <- element_names

        elements
}

#' Read a metadata file
#'
#' @param filename File name
#' @param id_col Identifier column in the file
#' @param sep File separator
#' @param stringsAsFactors Passed to \code{read.delim}
#'
#' @export
#'
#' @return output Data frame

read_metadata_yaml <- function(filename, id_col = NULL, sep = NULL, stringsAsFactors = FALSE) {
        if (is.null(sep)) {
                sep <- getSeparator(filename)
        }

        if (!file.exists(filename)) {
                stop(paste("Metadata file", filename, "does not exist."))
        } else {
                metadata <-
                        read.delim(
                                filename,
                                sep = sep,
                                check.names = FALSE,
                                header = TRUE,
                                stringsAsFactors = stringsAsFactors
                        )
        }

        if (!is.null(id_col)) {
                if (is.character(id_col) && !id_col %in% colnames(metadata)) {
                        stop(
                                paste0(
                                        "Metadata ID column (",
                                        id_col,
                                        ") does not exist in metadata ",
                                        paste(colnames(metadata), collapse = ","),
                                        " from file ",
                                        filename
                                )
                        )
                }

                metadata <- metadata[match(unique(metadata[[id_col]]), metadata[[id_col]]), ]
                rownames(metadata) <- metadata[[id_col]]
        }
        return(metadata)
}


#' Check one list is a subset of another and throw an error if not
#'
#' @param test_list  Test list
#' @param reference_list Reference list
#' @param test_list_name Name of test list for error
#' @param reference_list_name Name of reference list for error
#'
#' @return output Returns TRUE if check passes
#' @export

checkListIsSubset_yaml <- function(test_list,
                                                                                                                        reference_list,
                                                                                                                        test_list_name,
                                                                                                                        reference_list_name) {
        if (!all(test_list %in% reference_list)) {
                stop(
                        paste0(
                                "Not all ",
                                test_list_name,
                                " (",
                                paste(test_list, collapse = ","),
                                ") are available in the ",
                                reference_list_name,
                                " (",
                                paste(reference_list, collapse = ","),
                                ")"
                        )
                )
        }
        TRUE
}


#' Read and validate a contrasts file against sample metadata
#'
#' @param filename Contrasts file
#' @param samples Data frame of sample information
#' @param variable_column Column in contrasts file referencing sample sheet
#'   column
#' @param reference_column Column in contrast file referencing reference level
#'   of sample sheet variable
#' @param target_column Column in contrast file referencing target level of
#'   sample sheet variable
#' @param blocking_column Colon-separated column in contrast file referencing
#'   sample sheet variables to be used as blocking factors
#' @param convert_to_list Convert output to a list as used internally by
#'   shinyngs?
#'
#' @return output Validated contrasts data frame
#' @export

read_contrasts_yaml <- function(filename,
                                                                                                                                samples,
                                                                                                                                variable_column = "variable",
                                                                                                                                reference_column = "reference",
                                                                                                                                target_column = "target",
                                                                                                                                blocking_column = "blocking",
                                                                                                                                convert_to_list = FALSE) {

        # Read the contrasts depending on the file format (CSV or YAML)
        if (grepl("\\.csv$", filename)) {
                contrasts <- read_metadata_yaml(filename)
                contrast_cols <- c(variable_column, reference_column, target_column)
                if (!blocking_column %in% names(contrasts)) {
                        contrasts[[blocking_column]] <- NA
                }

                # Check contrast headers are as expected
                if (!all(contrast_cols %in% colnames(contrasts))) {
                        stop(paste("Contrasts file must contain all of", paste(contrast_cols, collapse = ", ")))
                }

        } else if (grepl("\\.ya?ml$", filename)) {
                contrasts_yaml <- yaml::read_yaml(filename)

                # Extract contrasts and structure into a data frame
                contrasts_list <- contrasts_yaml$contrasts
                contrasts <- data.frame(
                        id = sapply(contrasts_list, function(x) x$id),
                        variable = sapply(contrasts_list, function(x) x$comparison[1]),
                        reference = sapply(contrasts_list, function(x) x$comparison[2]),
                        target = sapply(contrasts_list, function(x) x$comparison[3]),
                        blocking = sapply(contrasts_list, function(x) ifelse(!is.null(x$blocking_factors), paste(x$blocking_factors, collapse = ";"), NA)),
                        stringsAsFactors = FALSE
                )

                # Check for missing fields
                if (any(is.na(contrasts$variable) | is.na(contrasts$reference) | is.na(contrasts$target))) {
                        stop("Contrasts file has missing values in key columns (variable, reference, target).")
                }
        } else {
                stop("Invalid file format. Please provide a CSV or YAML file.")
        }

        # Check contrast content is appropriate to sample sheet
        success <- checkListIsSubset_yaml(contrasts$variable, colnames(samples), "contrast variables", "sample metadata")

        # Check blocking variables, where supplied
        blocking <- unlist(lapply(contrasts[[blocking_column]], function(x) simpleSplit(x, ";")))
        blocking <- blocking[!is.na(blocking)]
        if (length(blocking) > 0) {
                success <- checkListIsSubset_yaml(blocking, colnames(samples), "blocking variables", "sample metadata")
        }

        # Ensure reference, target, and blocking values are valid for their variable
        for (i in 1:nrow(contrasts)) {
                var <- contrasts[i, variable_column]
                for (col in c(reference_column, target_column)) {
                        val <- contrasts[i, col]
                        if (is.na(val) || val == "") {
                                stop(paste("Missing value for", col, "in sample sheet"))
                        } else {
                                success <- checkListIsSubset_yaml(val, samples[[var]], "contrast levels", "sample metadata variable")
                        }
                }
        }

        # Convert contrasts to a list if requested
        if (convert_to_list) {
                contrasts <- apply(contrasts, 1, function(x) {
                        conlist <- split(unname(x), names(x))[names(x)]
                        rename <- c('variable' = 'Variable', 'reference' = 'Group.1', 'target' = 'Group.2')
                        rename_ind <- match(names(rename), names(conlist))
                        names(conlist)[rename_ind] <- rename
                        nonempty <- unlist(lapply(conlist, function(y) !(is.na(y) || is.null(y) || grepl("^\\s*$", y))))
                        conlist[nonempty]
                })
        }

        contrasts
}

        #' Call the various read/ validate methods for input data surrounding an experiment
#'
#' @param samples_metadata Sample metadata data frame
#' @param features_metadata Feature metadata data frame
#' @param assay_files List of assay matrices
#' @param contrasts_file Contrasts definition file
#' @param sample_id_col Column of sample metadata used for identifiers
#' @param feature_id_col Column of feature metadata used for identifiers
#' @param assay_names Optional comma-separated list of assay names
#' @param differential_results Optional list of differential stats files
#' @param pval_column P value column if differential stats files specified
#' @param qval_column Q value column if differential stats files specified
#' @param fc_column Fold change column if differential stats files specified
#' @param unlog_foldchanges Boolean- should fold changes in stats files be
#'   unlogged?
#'
#' @return output A named list with feature/ observation components
#' @export

validate_inputs_yaml <- function(samples_metadata,
                                                                                                                assay_files,
                                                                                                                contrasts_file = NULL,
                                                                                                                features_metadata = NULL,
                                                                                                                sample_id_col = "sample",
                                                                                                                assay_names = NULL,
                                                                                                                differential_results = NULL,
                                                                                                                feature_id_col = "gene_id",
                                                                                                                pval_column = "pval_column",
                                                                                                                qval_column = "qval_column",
                                                                                                                fc_column = "log2FoldChange",
                                                                                                                unlog_foldchanges = FALSE) {
        validated_parts <- list()

        # Read the sample (observation) - wise metadata

        print(paste(
                "Reading sample sheet at",
                samples_metadata,
                "with ID col",
                sample_id_col
        ))

        samples <- read_metadata_yaml(
                filename = samples_metadata,
                id_col = sample_id_col
        )
        validated_parts[[samples_metadata]] <- samples

        # Read feature-wise metadata if provided

        features <- NULL
        if (!is.null(features_metadata)) {
                print(paste(
                        "Reading feature metadata at",
                        features_metadata,
                        "with ID col",
                        feature_id_col
                ))

                features <- read_metadata_yaml(
                        filename = features_metadata,
                        id_col = feature_id_col
                )
                validated_parts[[features_metadata]] <- features
        }

        # Read the assay matrices

        assay_files <-
                stringsToNamedVector_yaml(
                        elements_string = assay_files,
                        simplify_files = FALSE,
                        prettify_names = FALSE
                )

        # Read the matrices while checking samples and features match columns and rows

        validated_parts[["assays"]] <- lapply(assay_files, function(x) {
                print(paste("Reading assay matrix", x, "and validating against samples and features (if supplied)"))

                mat <- read_matrix(
                        matrix_file = x,
                        sample_metadata = samples,
                        feature_metadata = features
                )
                print(paste("... ", x, "matrix good"))
                mat
        })

        # Read contrasts and check against sample info

        if (!is.null(contrasts_file)) {
                print("Reading contrast definitions and validating against sample sheet")
                validated_parts[[contrasts_file]] <- read_contrasts_yaml(contrasts_file, samples)
                print("... contrasts good")
        }

        if (!is.null(differential_results)) {
                contrast_stats_files <-
                        stringsToNamedVector_yaml(differential_results,
                                simplify_files = FALSE,
                                prettify_names = FALSE
                        )

                validated_parts[["differential_stats"]] <- lapply(contrast_stats_files, function(dsf) {
                        read_differential_yaml(
                                filename = dsf,
                                feature_id_column = feature_id_column,
                                pval_column = pval_column,
                                qval_column = qval_column,
                                fc_column = fc_column,
                                unlog_foldchanges = unlog_foldchanges
                        )
                })
        }

        validated_parts
}

        ########################################################



option_list <- list(
        make_option(
                c("-s", "--sample_metadata"),
                type = "character",
                default = NULL,
                help = "CSV-format sample metadata file."
        ),
        make_option(
                c("-i", "--sample_id_col"),
                type = "character",
                default = "sample",
                help = "Column in sample metadata used as sample identifier. Should be used to name columns of expression matrices, and duplicate rows will be removed based on this column."
        ),
        make_option(
                c("-f", "--feature_metadata"),
                type = "character",
                default = NULL,
                help = "TSV-format feature (often gene) metadata file."
        ),
        make_option(
                c("-j", "--feature_id_col"),
                type = "character",
                default = "gene_id",
                help = "Column in feature metadata used as feature identifier. Should be used to name columns of expression matrices."
        ),
        make_option(
                c("-e", "--assay_files"),
                type = "character",
                default = NULL,
                help = "Comma-separated list of TSV-format file expression matrix files."
        ),
        make_option(
                c("-c", "--contrasts_file"),
                type = "character",
                default = NULL,
                help = "CSV-format contrast file with variable,reference and target in the first 3 columns."
        ),
        make_option(
                c("-d", "--differential_results"),
                type = "character",
                default = NULL,
                help = "Tab-separated files containing at least fold change and p value, one for each row of the contrast file."
        ),
        make_option(
                c("-k", "--fold_change_column"),
                type = "character",
                default = "log2FoldChange",
                help = "Column in differential results files holding fold changes."
        ),
        make_option(
                c("-u", "--unlog_foldchanges"),
                action = "store_true",
                default = FALSE,
                help = "Set this option if fold changes should be unlogged."
        ),
        make_option(
                c("-p", "--pval_column"),
                type = "character",
                default = "padj",
                help = "Column in differential results files holding p values."
        ),
        make_option(
                c("-q", "--qval_column"),
                type = "character",
                default = "padj",
                help = "Column in differential results files holding q values/ adjusted p values."
        ),
        make_option(
                c("-o", "--output_directory"),
                type = "character",
                default = NULL,
                help = "Serialized R object which can be used to generate a shiny app."
        ),
        make_option(
                c("-t", "--separator"),
                type = "character",
                default = "\t",
                help = "Consistent separator for re-written files."
        )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check mandatory

mandatory <-
        c(
                "sample_metadata",
                "assay_files"
        )

missing_args <- mandatory[!mandatory %in% names(opt)]
if (length(missing_args) > 0) {
        stop(paste("Missing mandatory arguments:", paste(missing_args, collapse = ", ")))
}

library(shinyngs)

# validate_inputs_yaml() just wraps the parsing functions of shinyng, used by e.g.
# eselistfromConfig(). These functions are good for ensuring the consistency of
# FOM (feaure/ observation matrix) data.

validated_parts <- validate_inputs_yaml(
        samples_metadata = opt$sample_metadata,
        features_metadata = opt$feature_metadata,
        assay_files = opt$assay_files,
        assay_names = opt$assay_names,
        contrasts_file = opt$contrasts_file,
        sample_id_col = opt$sample_id_col,
        feature_id_col = opt$feature_id_col,
        differential_results = opt$differential_results,
        pval_column = opt$pval_column,
        qval_column = opt$qval_column,
        fc_column = opt$fold_change_column,
        unlog_foldchanges = opt$unlog_foldchanges
)

# If an output path is provided we can re-write the data, ensuring consistency
# of output formatting

if (! is.null(opt$output_directory)){

        dir.create(opt$output_directory, showWarnings = FALSE, recursive = TRUE)

        # Write the files back, but using the supplied separator

        write_table <- function(x, infile, suffix, na = 'NA'){
                file_basename <- tools::file_path_sans_ext(basename(infile))
                outfile <- file.path(opt$output_directory, paste(file_basename, suffix, 'tsv', sep = '.'))

                print(paste("...... writing", outfile))
                write.table(x, file = outfile, sep = opt$separator, quote = FALSE, row.names = FALSE, na = na)
        }

        # Write back the sample sheet, feature metadata and contrasts

        print("Writing basic data...")
        for (infile in c('sample_metadata', 'feature_metadata', 'contrasts_file')){
                filename <- opt[[infile]]
                if ((! is.null(filename)) && filename %in% names(validated_parts)){
                        write(paste("...", infile))

                        # Write contrasts file with empty strings for NAs in blocking
                        write_table(validated_parts[[filename]], filename, infile, na = ifelse(infile == 'contrasts_file', '', 'NA'))
                }
        }

        # Write back the matrices

        print("Writing matrices...")
        if ('assays' %in% names(validated_parts)){
                for (assay in names(validated_parts[['assays']])){
                        mat <- validated_parts[['assays']][[assay]]

                        # Add a column for row names
                        mat <- data.frame(feature_name = rownames(mat), mat, check.names = FALSE)
                        colnames(mat)[1] <- opt$feature_id_col

                        write_table(mat, assay, 'assay')
                }
        }

        # Write back the simplified differential results (if supplied)

        if ('differential_stats' %in% names(validated_parts)){
                for (ds in names(validated_parts[['differential_stats']])){
                        write_table(validated_parts[['differential_stats']][[ds]], ds)
                }
        }

}
