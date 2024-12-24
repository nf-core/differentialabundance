#!/usr/bin/env Rscript

## TODO: add make contrast arguments

# LOAD LIBRARIES ----------------------------------------------------
suppressWarnings(suppressMessages({
    library(tidyverse)
    library(optparse)
    library(yaml)
    library(jsonlite)
}))

# PARSE ARGUMENTS ---------------------------------------------------
## Generate list
option_list <- list(
    make_option(c("-y", "--yml"), type = "character", default = NULL,
                help = "Path to the models.yml file", metavar = "character"),
    make_option(c("-s", "--samplesheet"), type = "character", default = NULL,
                help = "Path to the samplesheet CSV file", metavar = "character"),
    make_option(c("-i", "--sample_id_col"), type = "character", default = "sample",
                help = "Column that contains sample identificators", metavar = "character")
)

## Parse command-line arguments
opt_parser <-
    OptionParser(
        option_list = option_list,
        description = "Validate a sample sheet against a YML file for model designs.",
        epilogue =
        "
        The process will look for variables and factors in the YML file and validate that they are present in the sample sheet.
        Will also look for unwanted characters within columns and colnames.
        Finally, will evaluate all model designs based on the formula and contrasts to find whether the models are full ranked or not.
        "
    )
opt <- parse_args(opt_parser)

## Validate input arguments
### Required arguments
if (is.null(opt$yml) || is.null(opt$samplesheet) ) {
    stop("'--yml' and '--samplesheet' arguments are required.", call. = FALSE)
}

### Collect parameters (easier for dev and testing)
path_yml         <- opt$yml             # "/workspace/differentialabundance/modules/local/validatemodel/tests/contrasts.yml" #
path_samplesheet <- opt$samplesheet     # "/workspace/differentialabundance/results/SRP254919.samplesheet_witherrors.csv" #
sample_column    <- opt$sample_id_col   # "sample" #

# LOAD FILES --------------------------------------------------------
## Load models.yml file
tryCatch({
    models <- read_yaml(path_yml)
    cat("Loaded YML file successfully.\n")
}, error = function(e) {
    stop("Error loading YML file: ", e$message)
})

## Load samplesheet CSV file
tryCatch({
    samplesheet <- readr::read_csv(path_samplesheet, show_col_types = FALSE)
    if ( !sample_column %in% colnames(samplesheet) ) {
        stop(paste0("Column ", sample_column, " not present in colnames(samplesheet)") )
    }
    cat("Loaded samplesheet successfully.\n")
}, error = function(e) {
    stop("Error loading samplesheet CSV: ", e$message)
})

# PARSE FACTORS/LEVELS FROM YML FILE --------------------------------
## Collect all variables and factors from the yml file
## * var_list:        collects all factors levels for each variable.
## * contrast_list:   collects info for each contrast

process_models <- function(models) {
    # Initialize lists to store the output
    contrasts_list <- list()
    var <- list()

    # Ensure models$models is not null or empty
    if (is.null(models$models) || length(models$models) == 0) {
        stop("models$models is null or empty. Please provide valid input.")
    }

    # Temporary storage for gathering contrasts per variable
    temp_var <- list()
    blocking_vars <- c()

    # Iterate through models (formula)
    for (FORMULA in models$models) {

        # Ensure FORMULA$contrasts is not null or empty
        if (is.null(FORMULA$contrasts) || length(FORMULA$contrasts) == 0) {
            stop("FORMULA$contrasts is null or empty. Unable to proceed with this model.")
        }

        # Iterate through contrasts for each model (formula)
        for (CONTRAST in FORMULA$contrasts) {

            # Validate CONTRAST$id
            if (is.null(CONTRAST$id)) {
                stop("Missing CONTRAST$id detected. Unable to proceed.")
            }

            if (CONTRAST$id %in% names(contrasts_list)) {
                stop(paste("Duplicate CONTRAST$id detected:", CONTRAST$id, ". Unable to proceed."))
            }

            # Get column name: first comparison's element
            name <- CONTRAST$comparison[1]

            # Check if name is valid
            if (is.null(name) || name == "") {
                stop("Invalid name detected in CONTRAST$comparison. Unable to proceed.")
            }

            # Get factors (rest of components)
            variables <- CONTRAST$comparison[2:length(CONTRAST$comparison)]

            # Populate contrasts_list for later model validation
            contrasts_list[[CONTRAST$id]] <- list(
                "formula" = FORMULA$formula,
                "variable" = name,
                "contrast" = variables,
                "blocking_factors" = CONTRAST$blocking_factors
            )

            # Gather contrasts per variable into temporary storage
            if (name %in% names(temp_var)) {
                temp_var[[name]] <- unique(c(temp_var[[name]], variables))
            } else {
                temp_var[[name]] <- variables
            }

            # Gather blocking factors
            if (!is.null(CONTRAST$blocking_factors)) {
                blocking_vars <- unique(c(blocking_vars, CONTRAST$blocking_factors))
            }
        }
    }

    # Consolidate gathered contrasts into `var`
    for (name in names(temp_var)) {
        if (name %in% names(var)) {
            var[[name]] <- unique(c(var[[name]], temp_var[[name]]))
        } else {
            var[[name]] <- temp_var[[name]]
        }
    }

    # Add blocking variables to `var`
    var[[ "blocking_factors" ]] <- blocking_vars

    # Return the populated lists
    return(list(
        contrasts_list = contrasts_list,
        var = var
    ))
}

## Process models
models_lists   <- process_models(models)
### Extract models contrasts
contrasts_list <- models_lists[["contrasts_list"]]
### Extract variables and levels
var            <- models_lists[["var"]]

## Print explicit messages
for (INDEX in 1:length(var)) {
    cat("\033[32mDetected '", names(var)[INDEX], "' variable with ", paste(var[[INDEX]], collapse = " "), " levels.\033[0m \n", sep = "")
}

# VALIDATE SAMPLESHEET BASED ON YML FILE ----------------------------
## Define function to compare YML and samplesheet
validate_model <- function(sample_column, variables, samplesheet) { # sample_column: string; variables: list; samplesheet: df

    #######################################################
    ##
    ## Function that takes a list of variables (and expected values) from the models.yml file that must be present in the sample sheet.
    ##
    ## The function evaluates that:
    ## * Every variable is present in the sample sheet.
    ## * Every variable has valid name
    ## * There are no NAs
    ## * All expected values are present
    ## * All samples contains a valid level.
    ##
    ## Finally, the function returns
    ## * An error vector containing all errors found.
    ## * A warnings vector
    ##
    ## If the error vector contains values, the script will end with a non-zero status.
    ##
    #######################################################

    ## Initialize errors vector
    errors <- c()

    ## Initialize warning errors
    warnings <- c()

    # Initialize vector to report continuos variables
    continuous <- c()

    ## Do not allow special characters
    undesired_chars <- "[^a-zA-Z0-9_.]"

    ## Check samplesheet names for invalid colnames
    df_colnames <- names(samplesheet)
    true_columns <- stringr::str_detect(df_colnames, pattern = regex(undesired_chars))

    if ( sum(true_columns) > 0 ) {
        errors <- c(errors,
            paste0("The following columns contain undesired characters: ", paste(df_colnames[true_columns], collapse = " "))
        )
    }

    ## Check that blocking variables exists and do not contain NAs, if they were specified
    blocking_factors <- c()
    if ( !is.null(variables$blocking_factors) ) {
        blocking_factors <- variables$blocking_factors

        for (VARIABLE in blocking_factors) {
            ## Check that the column exists
            if ( VARIABLE %in% colnames(samplesheet) ) {
                ## Check if there are NAs
                na_rows <- is.na( samplesheet[[ VARIABLE ]])
                if ( sum(na_rows) > 0 ) {
                    errors <- c(errors, paste0("Blocking factor ", VARIABLE, " contains NA/s in the following rows: ", paste(which(na_rows), collapse = " ")))
                }
            } else {
                errors <- c(
                    errors,
                    paste0("Blocking factor ", VARIABLE, " not present in sample sheet. Please check the compatibility between ", basename(path_yml), " and ", basename(path_samplesheet), " files."))
            }

            ## Alert about continuous variables
            if ( is.numeric(samplesheet[[ VARIABLE ]]) ) {
                continuous <- c(continuous, VARIABLE)
            } else {
                paste(VARIABLE, "not continuous")
            }
        }
    }

    ## Check variables and levels specified
    for (VARIABLE in names(variables)[ names(variables) != "blocking_factors" ] ) { ## Exclude blocking factors, they cannot be treated equally here

        ## Check that the column exists
        if ( VARIABLE %in% colnames(samplesheet) ) {

            ## Alert about continuous variables
            if ( is.numeric(samplesheet[[ VARIABLE ]]) ) {
                continuous <- unique(c(continuous, VARIABLE))
            } else {
                paste(VARIABLE, "is not continuous")
            }

            ## Check if there are NAs
            na_rows <- is.na( samplesheet[[ VARIABLE ]] )
            if ( sum(na_rows) > 0 ) {
                errors <- c(errors, paste0("Column ", VARIABLE, " contains NA/s in the following rows: ", paste(which(na_rows), collapse = " ")))
            }

            ## Check that the column data does not contain undesidered characters
            if ( sum( stringr::str_detect( samplesheet[[ VARIABLE ]], pattern = regex(undesired_chars) ), na.rm = TRUE) > 0 ) {
                errors <- c(errors,
                    paste0("Column ", VARIABLE, " contains undesired characters\n")
                )
            }

            ## Extract the expected levels according to the YML file
            expected_factors_levels <- c( rep(NA, length(variables[[ VARIABLE ]])) )

            ## Assign names
            names(expected_factors_levels) <- variables[[ VARIABLE ]]

            ## Extract the levels from the sample sheet
            obtained_factors_levels <- unique(samplesheet[[ VARIABLE ]])

            ## Check that all expected levels for this variable are effectively present
            for (LEVEL in names(expected_factors_levels) ) {
                expected_factors_levels[ LEVEL ] <- LEVEL %in% obtained_factors_levels
            }

            if( !all(expected_factors_levels) ) {
                errors <- c(
                    errors,
                    paste0(
                        "Missing factor levels for variable '", VARIABLE, "'. ",
                        "Present levels: ", paste( names( expected_factors_levels[ expected_factors_levels ]), collapse = " " ), ". ",
                        "Missing levels: ", paste( names( expected_factors_levels[!expected_factors_levels ]), collapse = " " )
                        )
                    )
            }

            ## Check whether there are levels in the table that are not present in the models definitions
            if (!all( obtained_factors_levels[!is.na(obtained_factors_levels)] %in% names(expected_factors_levels)) ) {
                warnings <- c(warnings,
                paste0( "The following labels are present in the samplesheet but are not part of the models definition and will be excluded: ", paste( obtained_factors_levels[!obtained_factors_levels %in% names(expected_factors_levels) ], collapse = " " ) ))
            }

        ## Report that a column is not even present in the table
        } else {
            errors <- c(
                errors,
                paste0("Column ", VARIABLE, " not present in sample sheet. Please check the compatibility between ", basename(path_yml), " and ", basename(path_samplesheet), " files."))
        }
    }

    ## Report if any variable was imported as continuous
    if ( !is.null(continuous)) {
        warnings <- c(warnings, paste0("The following continuous variables were detected or coerced into numeric: ", paste(continuous, collapse = "" )))
    }

    ## Report ERRORS and stop
    if (length(errors) > 0 ) {
        stop(cat("\033[1;31mSome errors where found while validating the samplesheet and models definitions:\n", paste(errors, collapse = "\n"), "\033[0m\n", sep = ""))
    }

    ## Generate validated phenotypic table
    tryCatch({
        ## Get desired columns, starting with the one containing sample identificators
        selected_columns <- c( sample_column, names(variables)[ names(variables) != "blocking_factors" ], blocking_factors)
        pheno_table      <- samplesheet %>% dplyr::select(all_of(selected_columns))

    }, error = function(e) {
        stop("Error generating validated phenotypic table: ", e$message)
    })

    ## Return list
    return(list( 'pheno_table' = pheno_table, 'warnings' = warnings) )

}

## Collect all outputs from function: errors must be empty, warnings can have messages
phenotable_warnings <- validate_model(sample_column, var, samplesheet)

## Cat warnings messages
cat("\033[1;33m", paste(phenotable_warnings[[ 2 ]], collapse = "\n"), "\033[0m\n")

## Get validated pheno table
pheno_table <- phenotable_warnings[[1]]

## FUNCTION TO CHECK THAT THE MODELS ARE FULL RANKED
check_model_contrasts <- function(contrasts_list, colData) {

    #######################################################
    ##
    ## Function that takes a list with models + variables, and a phenotypic table
    ##
    ## The function evaluates that:
    ## * The formulas are usable
    ## * The models generated are full ranked
    ##
    ## Finally, the function returns
    ## * An list with all designs
    ## * full ranked: true/false
    ##
    ## If the error vector contains values, the script will end with a non-zero status.
    ##
    #######################################################

    ### Initialize list for models descriptors
    design_list <- list()
    ### Iterate over each model
    for (model_name in names(contrasts_list)) {

        ## Extract model components
        model        <- contrasts_list[[model_name]]
        base_formula <- as.formula(model$formula)
        variable     <- model$variable
        contrast     <- model$contrast
        blocking     <- model$blocking_factor

        ## Check if blocking factor is already in the formula
        if (!is.null(blocking)) {
            ## Get terms from the formula
            formula_terms <- all.vars(base_formula)

            ## If blocking factor is not already in the formula, add it
            if (!(blocking %in% formula_terms)) {
                updated_formula <- as.formula(paste(deparse(base_formula), "+", blocking))
            } else {
                updated_formula <- base_formula
            }
        } else {
            updated_formula <- base_formula
        }

        cat("\nChecking:", model_name, "\nFormula:", deparse(updated_formula), "\n")

        # Build the design matrix
        design_matrix <- model.matrix(updated_formula, data = colData)

        # Check the rank of the design matrix
        rank <- qr(design_matrix)$rank
        expected_rank <- ncol(design_matrix)

        ## Add results to list
        design_list[[ length(design_list) + 1 ]] <- list(
            "formula"       = deparse(updated_formula),
            "design_matrix" = design_matrix,
            "full_rank"     = rank == expected_rank
        )
        names(design_list)[ length(design_list) ]  <- model_name
    }

    for (DESIGN in names(design_list)) {
        cat("Design ", DESIGN,
            if( design_list[[ DESIGN ]]$full_rank ) { " is " } else { " is not "},
            "full ranked\n", sep = ""
        )
    }

    return(design_list)
}

design_results <- check_model_contrasts(contrasts_list, pheno_table)

# EXPORT DATA -------------------------------------------------------
## Export pheno table
write_csv(pheno_table, "pheno_table.csv")

## Export designs to JSON
jsonlite::write_json(design_results, path = "designs.json", pretty = TRUE, auto_unbox = TRUE)

## Export warnings
if ( !is.null(phenotable_warnings[[ 2 ]]) ) {
    jsonlite::write_json(phenotable_warnings[[ 2 ]], path = "warnings.json", pretty = TRUE)
}

## Export RData
save.image("models.RData")
