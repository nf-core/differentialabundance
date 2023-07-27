#!/usr/bin/env Rscript

#Libraries

library(qusage)
library(GSVA)
library(org.Hs.eg.db)
library(magrittr)
library(msigdbr)

####Commandline Argument parsing###
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("Usage: compute_GSVA.R <scores.yaml> <counts.tsv> species", call.=FALSE)
}
####Uncomment if debugging

gene_sets_to_compute <- args[1]
counts_input <- args[2]
species <- args[3]
output_dir <- getwd()
counts <- read.table(counts_input, sep="\t", check.names=F, header=T)
output_basename <- paste0(basename(tools::file_path_sans_ext(counts_input)), "_gsva_results.tsv")

hallmark_gene_sets <- msigdbr::msigdbr(
  species = species, # Can change this to what species you need
  category = gene_sets_to_compute # Only hallmark gene sets
)
hallmarks_list <- split(
  hallmark_gene_sets$entrez_gene, # The genes we want split into pathways
  hallmark_gene_sets$gs_name # The pathways made as the higher levels of the list
)

mapped_df <- data.frame(
  "entrez_id" = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Hs.eg.db,
    keys = counts$gene_id,
    # Replace with the type of gene identifiers in your data
    keytype = "ENSEMBL",
    # Replace with the type of gene identifiers you would like to map to
    column = "ENTREZID",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a Entrez gene identifier,
  # drop that from the data frame
  dplyr::filter(!is.na(entrez_id)) %>%
  # Make an `Ensembl` column to store the row names
  tibble::rownames_to_column("Ensembl") %>%
  # Now let's join the rest of the expression data
  dplyr::inner_join(counts, by = c("Ensembl" = "gene_id"))

gene_means <- rowMeans(mapped_df %>% dplyr::select(-Ensembl, -entrez_id))

# Let's add this as a column in our `mapped_df`.
mapped_df <- mapped_df %>%
  # Add gene_means as a column called gene_means
  dplyr::mutate(gene_means) %>%
  # Reorder the columns so `gene_means` column is upfront
  dplyr::select(Ensembl, entrez_id, gene_means, dplyr::everything())

filtered_mapped_df <- mapped_df %>%
  # Sort so that the highest mean expression values are at the top
  dplyr::arrange(dplyr::desc(gene_means)) %>%
  # Filter out the duplicated rows using `dplyr::distinct()`
  dplyr::distinct(entrez_id, .keep_all = TRUE)

filtered_mapped_matrix <- filtered_mapped_df %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-Ensembl, -gene_means) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("entrez_id") %>%
  # Now we can convert our object into a matrix
  as.matrix()
gsva_results <- gsva(
  filtered_mapped_matrix,
  hallmarks_list,
  method = "gsva",
  kcdf = "Gaussian",
  # Minimum gene set size
  min.sz = 15,
  # Maximum gene set size
  max.sz = 500,
  # Compute Gaussian-distributed scores
  mx.diff = TRUE,
  # Don't print out the progress bar
  verbose = TRUE
)

gsva_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pathway") %>%
  readr::write_tsv(file.path(
    output_dir,
    output_basename
  ))
