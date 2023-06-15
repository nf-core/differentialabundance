#!/usr/bin/env Rscript

############################
# FUNCTIONS
############################
# From affy/justRMA (pinin4fjords)
#  Parse out options from a string without recourse to optparse
# @param x Long-form argument list like --opt1 val1 --opt2 val2
# return named list of options and values similar to optparse
parse_args <- function(x){
  args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
  args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

  # Ensure the option vectors are length 2 (key/ value) to catch empty ones
  args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

  parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
  parsed_args[! is.na(parsed_args)]
}
# From affy/justRMA (pinin4fjords)
# Round numeric dataframe columns to fixed decimal places by applying
# formatting and converting back to numerics
# @param dataframe A data frame
# @param columns Which columns to round (assumes all of them by default)
# @param digits How many decimal places to round to?
# @return output Data frame
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

############################
# PARSE PARAMS FROM NEXTFLOW
############################

opt <- list(
  querygse = '$querygse'
)
args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
  if (! ao %in% names(opt)){
    stop(paste("Invalid option:", ao))
  }else{
    opt[[ao]] <- args_opt[[ao]]
  }
}

############################
# MAIN
############################

library(GEOquery)

# fetch data for GSE number
eset <- getGEO(opt\$querygse)[[1]]

# write probeset annotation
write.table(fData(eset)[,c('ID','Entrez_Gene_ID','Symbol','Definition')],
            paste0(opt\$querygse,'.annotation.tsv'),
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)


# if data is not log scale, transform it as needed for limma downstream
if(max(exprs(eset)) > 20) { # a bit dirty, needs proper solution later...
  exprs(eset)[exprs(eset) <= 0] <- .001
  exprs(eset) <- log2(exprs(eset))
}

output_prefix <- '$task.ext.prefix'
saveRDS(eset, file = paste0(output_prefix, 'eset.rds'))

# write intensity matrix (normalised)
write.table(
  data.frame(
    probe_id = rownames(eset),
    round_dataframe_columns(as.data.frame(exprs(eset))),
    check.names = FALSE
  ),
  file = paste0(output_prefix, 'matrix.tsv'),
  col.names = TRUE, row.names = FALSE,
  sep = '\t', quote = FALSE
)


############################
# LOG SESSION AND VERSIONS
############################


sink("R_sessionInfo.log")
print(sessionInfo())
sink()

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
geoquery.version <- as.character(packageVersion("GEOquery"))

writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r.version),
    paste('    bioconductor-:', geoquery.version)
  ),
  'versions.yml')
