# nf-core/differentialabundance: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/differentialabundance/usage](https://nf-co.re/differentialabundance/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

Differential analysis is a common task in a variety of use cases. In essence, all these use cases entail taking an input matrix containing features (e.g. genes) and observations (e.g. samples), and comparing groups of observations in all or a subset of the features. The feature/ observation language here reflects our hope that this workflow will extend in future to encompass a variety of applications where an assumption of gene vs sample may not be a valid one - though that is the application to which the first release will apply.

With the above in mind, running this workflow requires:

- a set of abundance values. This can be:
  - (for RNA-seq or MaxQuant proteomics measurements): a matrix of quantifications with observations by column and features by row
  - (for Affymetrix microarrays): a tar'd archive of CEL files
- a description of the observations such as a sample sheet from RNA-seq analysis
- a description of the features, for our initial RNA-seq application this can be simply the GTF file from which gene annotations can be derived. For Affymetrix arrays this can be derived from the array platform annotation package automatically. Skip for MaxQuant. You can also supply your own table.
- a specification of how the matrix should be split, and how the resulting groups should be compared

## Observations (samplesheet) input

```bash
--input '[path to samplesheet file].(csv|tsv)'
```

The samplesheet file can be tab or comma separated. This may well be the same sample sheet used to generate the input matrix. For example, in RNA-seq this might be the same sample sheet, perhaps derived from [fetchngs](https://github.com/nf-core/fetchngs), that was input to the [RNA-seq workflow](https://github.com/nf-core/rnaseq). It may be necessary to add columns that describe the groups you want to compare. The columns that the pipeline requires are:

- a column listing the sample IDs (must be the same IDs as in the abundance matrix), in the example below it is called `sample`. For some study_types, this column might need to be filled in with file names, e.g. when doing an affymetrix analysis.
- one or more columns describing conditions for the differential analysis. In the example below it is called `condition`
- optionally one or more columns describing sample batches or similar which you want to be considered in the analysis. In the example below it is called `batch`

For example:

```csv
sample,fastq_1,fastq_2,condition,replicate,batch
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,control,1,A
CONTROL_REP2,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,control,2,B
CONTROL_REP3,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,control,3,A
TREATED_REP1,AEG588A2_S1_L002_R1_001.fastq.gz,AEG588A2_S1_L002_R2_001.fastq.gz,treated,1,B
TREATED_REP2,AEG588A2_S1_L003_R1_001.fastq.gz,AEG588A2_S1_L003_R2_001.fastq.gz,treated,2,A
TREATED_REP3,AEG588A2_S1_L004_R1_001.fastq.gz,AEG588A2_S1_L004_R2_001.fastq.gz,treated,3,B
```

### Affymetrix arrays

Abundances for Affy arrays are provided in CEL files within an archive. When creating sample sheets for Affy arrays, it's crucial to include a column that specifies which file corresponds to each sample. This file column is essential for linking each sample to its corresponding data file, as shown in the example below:

```
"file","id","name","patient","phenotype"
"GSM1229341_Gudjohnsson_001_6690_PP.CEL.gz","GSM1229341","p6690_PP","6690","lesional"
"GSM1229342_Gudjohnsson_002_6690_PN.CEL.gz","GSM1229342","p6690_PN","6690","uninvolved"
"GSM1229343_Gudjohnsson_003_7450_PN.CEL.gz","GSM1229343","p7450_PN","7450","uninvolved"
"GSM1229344_Gudjohnsson_004_7450_PP.CEL.gz","GSM1229344","p7450_PP","7450","lesional"
"GSM1229345_Gudjohnsson_005_7912_PP.CEL.gz","GSM1229345","p7912_PP","7912","lesional"
"GSM1229346_Gudjohnsson_006_7912_PN.CEL.gz","GSM1229346","p7912_PN","7912","uninvolved"
"GSM1229347_Gudjohnsson_007_8470_PP.CEL.gz","GSM1229347","p8470_PP","6690","lesional"
"GSM1229348_Gudjohnsson_008_8470_PN.CEL.gz","GSM1229348","p8470_PN","6690","uninvolved"
```

The `file` column in this example is used to specify the data file associated with each sample, which is essential for data analysis and interpretation.

## Abundance values

### RNA-seq and similar

```bash
--matrix '[path to matrix file].(csv|tsv)'
```

This is a numeric matrix file, comma or tab-separated, with features as rows and observations in columns. The features correspond to the supplied feature set. The parameters `--observations_id_col` and `--features_id_col` define which of the associated fields should be matched in those inputs.

#### Outputs from nf-core/rnaseq and other tximport-processed results

The nf-core RNAseq workflow incorporates [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html) for producing quantification matrices. From [version 3.12.2](https://github.com/nf-core/rnaseq/releases/tag/3.13.2), it additionally provides transcript/gene length matrices which can be directly consumed by DESeq2 to model length bias across samples.

To use this approach, include the corresponding lengths file with the **raw counts**:

```bash
--matrix 'salmon.merged.gene_counts.tsv' \
--feature_length_matrix 'salmon.merged.gene_lengths.tsv'
```

Without the feature lengths, for instance in earlier rnaseq workflow versions, follow the second recommendation in the [tximport documentation](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#Downstream_DGE_in_Bioconductor):

> "Use the tximport argument `countsFromAbundance='lengthScaledTPM'` or `'scaledTPM'`, then employ the gene-level count matrix `txi$counts` directly in downstream software, a method we call 'bias corrected counts without an offset'"

This aligns with the **gene_counts_length_scaled.tsv** or **gene_counts_scaled.tsv** matrices in the rnaseq workflow.

It is important to note that the documentation advises:

> "Do not manually pass the original gene-level counts to downstream methods without an offset."

So we **do not recommend** raw counts files such as `salmon.merged.gene_counts.tsv` as input for this workflow **except** where the transcript/gene lengths are also provided.

### MaxQuant intensities

```bash
--matrix '[path to matrix file]'
```

This is the proteinGroups.txt file produced by MaxQuant. It is a tab-separated matrix file with a column for every observation (plus additional columns for other types of measurements and information); each row contains these data for a set of proteins. The parameters `--observations_id_col` and `--features_id_col` define which of the associated fields should be matched in those inputs. The parameter `--proteus_measurecol_prefix` defines which prefix is used to extract those matrix columns which contain the measurements to be used. For example, the default `LFQ intensity ` will indicate that columns like LFQ intensity S1, LFQ intensity S2, LFQ intensity S3 etc. are used (one whitespace is automatically added if necessary).

### Generic pre-scaled matrix

```bash
-profile generic_matrix
--matrix '[path to matrix file].(csv|tsv)'
```

Use this profile when you already have an abundance matrix on the appropriate scale (e.g. log-transformed microarray intensities, log-CPM values, or any other pre-normalised data) and do not need RNA-seq-specific preprocessing such as DESeq2's VST/rlog or voom. The matrix is passed directly to Limma for differential analysis without any additional transformation.

As with RNA-seq input, the matrix should have features as rows and observations as columns. Set `--observations_id_col` and `--features_id_col` to match the relevant identifier columns in the samplesheet and matrix, respectively. Feature annotations can be provided via `--features` or, if omitted, the workflow will fall back to using the matrix row identifiers directly.

For mixed-effects models, use `-profile generic_matrix_dream` instead, which runs DREAM rather than Limma.

### Affymetrix microarrays

```bash
--affy_cel_files_archive '[path to an archive of CEL files]'
```

This is an archive of CEL files as frequently found in GEO.

### Use SOFT matrices

Alternatively, the user may want to work with SOFT matrices. In this case, setting

`--study_type geo_soft_file` and `--querygse [GSE study ID]`

enables the pipeline to download normalised SOFT matrices automatically (note that even though Affymetrix arrays are also supported in the SOFT matrix track, it is recommended to work from CEL files in this case).

As for other platforms You may subset the metadata features used in reporting etc. e.g. for GPL570 (Affymetrix Plus 2.0 arrays) this could be done with

```
--features_metadata_cols ID,Entrez_Gene_ID,Symbol,Definition
```

Full list of features metadata are available on GEO platform pages.

## Contrasts file

The contrasts file references the observations file to define groups of samples to compare. It can be provided in **either** CSV/TSV or YAML format using a single `--contrasts` parameter. The file format is inferred from the extension (`.csv`, `.tsv`, `.yml`, `.yaml`).

### CSV/TSV contrasts file

```bash
--contrasts '[path to contrasts file].(csv|tsv)'
```

Based on the sample sheet above we could define contrasts as indicated below:

```csv
id,variable,reference,target,blocking
condition_control_treated,condition,control,treated,
condition_control_treated_blockrep,condition,control,treated,replicate;batch
```

The necessary fields in order are:

- `id` - an arbitrary identifier, will be used to name contrast-wise output files
- `variable` - which column from the observations information will be used to define groups
- `reference` - the base/ reference level for the comparison. If features have higher values in this group than target they will generate negative fold changes
- `target` - the target/ non-reference level for the comparison. If features have higher values in this group than the reference they will generate positive fold changes

You can optionally supply:

- `blocking` - semicolon-delimited, any additional variables (also observation columns) that should be modelled alongside the contrast variable
- `exclude_samples_col` and `exclude_samples_values` - the former being a valid column in the samples sheet, the latter a semicolon-delimited list of values in that column which should be used to select samples prior to differential modelling. This is helpful where certain samples need to be excluded prior to analysis of a given contrast.

### YAML contrasts file format

```bash
--contrasts '[path to YAML contrasts file]'
```

Based on the sample sheet above we could define YAML contrasts like:

```yaml
contrasts:
  - id: condition_control_treated
    comparison: ["condition", "control", "treated"]
  - id: condition_control_treated_blockrep
    comparison: ["condition", "control", "treated"]
    blocking_factors: ["replicate"]
```

The necessary fields in order are:

- `id` - an arbitrary identifier, will be used to name contrast-wise output files
- `comparison`(respectively):
- - `variable` - which column from the observations information will be used to define groups
- - `reference` - the base/ reference level for the comparison. If features have higher values in this group than target they will generate negative fold changes
- - `target` - the target/ non-reference level for the comparison. If features have higher values in this group than the reference they will generate positive fold changes
- `blocking_factors` - Any additional variables (also observation columns) that should be modelled alongside the contrast variable
- `exclude_samples_col` and `exclude_samples_values` - the former being a valid column in the samples sheet, the latter a list of values in that column which should be used to select samples prior to differential modelling. This is helpful where certain samples need to be excluded prior to analysis of a given contrast.

Additionally, the YAML contrasts also supports formula based model definitions:

```yaml
contrasts:
  - id: condition_control_treated
    formula: "~ condition"
    make_contrasts_str: "conditiontreated"
  - id: condition_control_treated_blockrep
    formula: "~ condition + replicate"
    make_contrasts_str: "conditiontreated"
```

The necessary fields in order are:

- `formula` - A string representation of the model formula. It is used to build the design matrix.
- `make_contrasts_str` - An explicit literal contrast string (e.g., "treatmenthND6 - treatmentmCherry") that is passed directly to [`limma::makeContrasts()`](https://rdrr.io/bioc/limma/man/makeContrasts.html) in `VARIANCEPARTITION_DREAM`, `LIMMA_DIFFERENTIAL` and `DESEQ2_DIFFERENTIAL`. The parameter names must be syntactically valid variable names in R (see [`make.names`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/make.names.html)). This field provides full control for complex designs. Requires `formula`.

> [!IMPORTANT]
>
> - YAML contrast definitions using `comparison` **and** those using `formula` are both supported by all differential methods.
> - They can be freely mixed within the same `contrasts:` list in a single YAML file (i.e. some contrasts may use `comparison` while others use `formula`).

> [!NOTE]
>
> #### Notes on `make_contrasts_str`
>
> This string must match exactly the name of the coefficient in the model matrix as generated by the specified `formula`. It is passed to `limma::makeContrasts()` without modification. For example:
>
> - `formula: "~ condition"` will generate model coefficients like `conditiontreated` (if `control` is the reference).
> - Then, `make_contrasts_str: "conditiontreated"` selects that coefficient for testing.
>
> This gives full control over the contrast definition but requires understanding of the model matrix.
> Some downstream applications (e.g. `GSEA_GSEA`, `SHINYNGS_APP`) do not support formula-based contrasts as they require a `meta.variable`.

Beyond the basic one-factor comparison, the YAML contrasts format supports advanced experimental designs through the use of interaction terms and custom contrast strings. These are particularly useful in multifactorial experiments where the effect of one variable may depend on the level of another (e.g. genotype × treatment). To model an interaction between genotype and treatment, use a formula like `~ genotype * treatment`, which expands the yaml to:

```yaml
contrasts:
  - id: genotype_WT_KO_treatment_Control_Treated
    formula: "~ genotype * treatment"
    make_contrasts_str: "genotypeKO.treatmentTreated"
```

To facilitate constructing and validating such models and contrast strings, consider using the [`ExploreModelMatrix`](https://www.bioconductor.org/packages/release/bioc/html/ExploreModelMatrix.html) Shiny app to have visual inspection of the design matrix and interactive contrast building. Another helpful resource is the [guide to creating design matrices for gene expression experiments](https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html).

## Feature annotations

### GTF file

```bash
--gtf '[path to gtf file]'
```

This is usually the easiest way to supply annotations for RNA-seq features. It should match the GTF used in nf-core/rnaseq if that workflow was used to produce the input expression matrix. Skip for MaxQuant.

### Annotation package identifiers for Affymetrix arrays

For `-profile affy`, default behaviour is to derive an annotation table while running the affy/justrma module based on the CDF name discovered there.

### Your own features, or no features

To override the above options, you may also supply your own features table as a TSV:

```bash
--features '[path to features TSV]'
```

By default, if you don't provide features, for non-array data the workflow will fall back to attempting to use the matrix itself as a source of feature annotations. For this to work you must make sure to set the `features_id_col`, `features_name_col` and `features_metadata_cols` parameters to the appropriate values, for example by setting them to 'gene_id' if that is the identifier column on the matrix. This will cause the gene ID to be used everywhere rather than more accessible gene symbols (as can be derived from the GTF), but the workflow should run. Please use this option for MaxQuant analysis, i.e. do not provide features.

## Reproducibility and random seeds

Some supported methods have stochastic components. You can now set a pipeline-wide seed with:

```bash
--seed 1234
```

When provided, this value is passed to supported stochastic methods, currently including DESeq2, limma, DREAM and GSEA, so reruns can be reproduced consistently across the pipeline.

If `--seed` is left unset (the default), the pipeline does not set a seed for these methods, so stochastic steps may remain non-deterministic.

## Mixed-effects models with DREAM

When your design has repeated measurements on the same biological unit (e.g. multiple samples per donor, technical replicates within subjects, time-course on the same individuals), the assumption of independent observations that DESeq2 and standard Limma rely on is violated. The pipeline supports DREAM (from the `variancePartition` R package), which extends Limma's linear-modelling framework with random effects so that within-subject correlation can be modelled explicitly.

To run DREAM, select one of the dedicated profiles:

- `-profile rnaseq_dream` for raw RNA-seq counts (voom transformation applied within the module).
- `-profile generic_matrix_dream` for a pre-scaled abundance matrix.

These profiles set `differential_method = 'dream'` and apply DREAM-specific defaults; do not override `--differential_method` on the command line.

Random effects are specified inside the [YAML contrasts file](#yaml-contrasts-file-format) using the `formula` field, with `(1 | variable)` syntax following `lme4` conventions. For example, with a `donor` column in the samplesheet:

```yaml
contrasts:
  - id: condition_control_treated
    formula: "~ condition + (1 | donor)"
    make_contrasts_str: "conditiontreated"
```

DREAM-specific parameters (`--dream_p_value`, `--dream_lfc`, `--dream_ddf`, `--dream_reml`, `--dream_apply_voom`, `--dream_adjust_method`, and others) control p-value adjustment, effect-size thresholds, the degrees-of-freedom approximation, REML vs ML fitting, and whether voom is applied. See the [parameters page](https://nf-co.re/differentialabundance/parameters) (the "Differential dream specific options" section) for the full list and defaults.

For background on the method, see the [DREAM paper (Hoffman & Roussos 2021)](https://pubmed.ncbi.nlm.nih.gov/32730587/) and the upstream [variancePartition documentation](https://bioconductor.org/packages/release/bioc/html/variancePartition.html).

## Analysis modes

The pipeline supports two modes of operation, each with well-defined parameter precedence rules:

> [!NOTE]
> The `study_type` parameter describes the input data format category only (for example matrix input, Affymetrix CEL archives, MaxQuant tables, or GEO SOFT retrieval). It is used for input validation and routing of format-specific preprocessing steps. It does not configure which statistical methods are run.

Method selection is configured through analysis profiles (for example `-profile rnaseq`, `-profile rnaseq_limma`, `-profile rnaseq_dream`) or explicit method parameters. Because of this separation, an RNA-seq profile and `study_type: rnaseq` are related defaults but conceptually different settings.

For matrix-style tabular data that is not RNA-seq, you can use `study_type: generic_matrix` (equivalent to `study_type: rnaseq` in pipeline behaviour). The `rnaseq` name is retained for backwards compatibility.

### 1. Single-Run Mode (Config Profiles) — recommended for production

For standard production use, select an **analysis profile** that bundles the correct study type, differential method, and output settings. CLI flags override profile parameters, following standard Nextflow precedence:

```
CLI flags / -params-file > profile > defaults
```

For example, to run an RNA-seq analysis with DESeq2:

```bash
nextflow run nf-core/differentialabundance \
    -profile rnaseq,docker \
    --input samplesheet.csv \
    --contrasts contrasts.yaml \
    --matrix counts.tsv \
    --gtf genes.gtf \
    --outdir results/
```

You can override any profile parameter from the command line. For example, to switch from DESeq2's default variance-stabilising transform to rlog:

```bash
nextflow run nf-core/differentialabundance \
    -profile rnaseq,docker \
    --deseq2_vs_method rlog \
    --input samplesheet.csv \
    --contrasts contrasts.yaml \
    --matrix counts.tsv \
    --gtf genes.gtf \
    --outdir results/
```

> [!WARNING]
> Do not override `--differential_method` on the command line when using an analysis profile. Each profile also sets method-specific parameters (e.g. logFC and p-value column names). To change methods, use the appropriate profile (e.g. `-profile rnaseq_limma`).

See the [Profiles](#-profile) section for the full list of available analysis profiles.

### 2. Multi-Run Mode (Paramsheet) — for exploration

For running multiple analysis configurations in parallel (e.g. comparing DESeq2 vs Limma, or testing different enrichment tools), provide a custom **paramsheet** — a YAML file where each entry defines a set of parameter overrides:

```bash
nextflow run nf-core/differentialabundance \
    -profile docker \
    --paramsheet my_configs.yaml \
    --input samplesheet.csv \
    --contrasts contrasts.yaml \
    --matrix counts.tsv \
    --gtf genes.gtf \
    --outdir results/
```

> [!WARNING]
> In multi-run mode, paramsheet parameters take precedence over CLI flags. The priority order is: `paramsheet > CLI flags > defaults`. This is by design — each paramsheet entry fully defines its configuration. CLI parameters are still used for values not specified in the paramsheet (e.g. `--input`, `--outdir`).

To run only a subset of configurations from your paramsheet, use `--paramset_name` with a comma-separated list:

```bash
--paramset_name deseq2_rnaseq,limma_rnaseq
```

> [!NOTE]
> The pipeline does **not** ship a default paramsheet. Users provide their own when using multi-run mode.

A paramsheet entry looks like this:

```yaml
- paramset_name: my_deseq2_run
  study_type: rnaseq
  differential_method: deseq2
  # ... additional parameter overrides
```

Each entry must include a unique `paramset_name`. Entries can override any pipeline parameter.

## Working with the output Quarto file

The pipeline produces a Quarto document file which, if you're proficient in R, you can use to tweak the report after it's generated.

> [!NOTE]
>
> If you need the same customisations repeatedly we would recommend you supply your own templates using the `report_file` parameter. Multiple templates can be supplied as a comma separated list.

To work with Quarto document files you will need an R environment with the required packages installed, such as the ShinyNGS R module [installed](https://github.com/pinin4fjords/shinyngs#installation), since it supplies a lot of the accessory plotting functions etc that you will need. An editor such as VS Code is recommended for viewing and editing the files. The way you will do this may depend on your exact systems, but for example:

### 1. Create a conda environment with Shinyngs and activate it

```bash
conda create -f ./assets/report_environment.yml -n report_environment
conda activate report_environment
```

### 2. Unzip the report archive

Now, unzip the report archive, and change the directory to that location.

### 3. Render Quarto report

Once the environment is active, you should have everything required to execute the code chunks and render the HTML report.

To render the report and open a live preview in your browser:

```bash
quarto preview differentialabundance_report.qmd --execute-params params.yml
```

For a better understanding of Quarto, you can go to the [comprehensive guide](https://quarto.org/docs/guide/).

## Shiny app generation

The pipeline is capable of building, and even deploying (to [shinyapps.io](https://www.shinyapps.io/)) for you a Shiny app built with [ShinyNGS](https://github.com/pinin4fjords/shinyngs). There is a basic example running [here](https://pinin4fjords.shinyapps.io/tester/) which shows what this might look like.

This is enabled with:

```bash
--shinyngs_build_app true
```

... which is the default. By default the app is not deployed, but just output to the output folder under `shinyngs_app/[study_name]`.

You have 3 choices in running that application:

1. Run locally
2. Have shinyapps.io host it for you
3. Host on a Shiny server

### 1. Run locally

You can start the application locally (in an environment where [ShinyNGS](https://github.com/pinin4fjords/shinyngs) is installed) like:

```bash
cd [output directory]/[study id]
Rscript app.R
```

This will give you a local URI to access in your browser:

```
Listening on http://127.0.0.1:3326
```

### 2. Shinyapps.io deployment

shinyapps.io is a hosting solution supplied by Posit (formerly RStudio) which gives you quick and easy access to hosting for Shiny applications. There is a free tier, though you'll have to pay for features such as authentication and improved resources.

You can upload your app to shinyapps.io youself, or deploy directly to shinyapps.io with this workflow, for which a few things need to happen:

#### Account and app setup

At https://www.shinyapps.io/, create an account, add a token (via Account -> Tokens) and note your secret and token.

You let Nextflow know about these via secrets:

```bash
nextflow secrets set SHINYAPPS_TOKEN [token]
nextflow secrets set SHINYAPPS_SECRET [secret]
```

#### Configuration

You then need to activate the deployment in your parameters, and supply both your account name and an app name:

```bash
--shinyngs_deploy_to_shinyapps_io \
--shinyngs_shinyapps_account '[account name]' \
--shinyngs_shinyapps_app_name '[app name]'
```

With this configuration in place deployment should happen automatically every time you run your workflow.

### 3. Run your own Shiny server

There is also a [Shiny server application](https://posit.co/download/shiny-server/), which you can install on your own infrastruture and use to host applications yourself.

## Gene set enrichment analysis

Currently, three tools can be used to do gene set enrichment analysis.

### GSEA

[GSEA](https://www.gsea-msigdb.org/gsea/index.jsp) tests for differential genes from within a user-provided set of genes; this requires a GMT or GMX file. The following example shows how to enable this:

```bash
--functional_method gsea \
--gene_sets_files gene_sets.gmt
```

### gProfiler2

The [gprofiler2](https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html) package can be used to test which pathways are enriched in the sets of differential genes produced by the the DESeq2 or limma modules. It is an R interface for the gprofiler webtool. In the simplest form, this feature can be enabled with the parameters from the following example:

```bash
--functional_method gprofiler2 \
--gprofiler2_organism mmusculus
```

If gene sets have been specified to the workflow via `--gene_sets_files` these are used by default. Specifying `--gprofiler2_organism` (mmusculus for Mus musculus, hsapiens for Homo sapiens etc.) will override those gene sets with gprofiler's own for the relevant species. `--gprofiler2_token` will override both options and use gene sets from a previous gprofiler run.

By default the analysis will be run with a background list of genes that passed the abundance filter (i.e. those genes that actually had some expression); see for example ["Multiple sources of bias confound functional enrichment analysis of global -omics data"](https://doi.org/10.1186/s13059-015-0761-7) for why this is advisable. You can provide your own background list with `--gprofiler2_background_file background.txt`or if you want to not use any background, set `--gprofiler2_background_file false`.

Check the [pipeline webpage](https://nf-co.re/differentialabundance/parameters#gprofiler2) for a full listing of the relevant parameters.

### Decoupler

[Decoupler](https://decoupler-py.readthedocs.io/en/latest/index.html) is a Python function that infers biological regulator activities—such as transcription factor or pathway activity—from omics data using multiple statistical enrichment methods. It takes as input a gene expression matrix and a prior knowledge network linking regulators to target genes, and applies one or more methods (e.g., ULM, MLM, wsum) to estimate regulator activity scores across samples. The function supports optional consensus scoring and outputs method-specific activity estimates and p-values, making it a versatile tool for activity inference in both bulk and single-cell datasets. If you want to see the full list of available methods and functions, refer to the function's [official guide](https://decoupler.readthedocs.io/en/v1.9.2/generated/decoupler.decouple.html#decoupler.decouple).

This tool is turned off by default, the following example shows how to enable it:

```bash
--functional_method decoupler \
--decoupler_network network.tsv
```

#### Input Files

Decoupler needs a matrix (mat) of molecular readouts (gene expression, logFC, p-values, etc.) and a network (net) that relates target features (genes, proteins, etc.) to “source” biological entities (pathways, transcription factors, molecular processes, etc.).

- The matrix will be taken from the results of the differential expression analysis performed by DESeq2, Limma, or DREAM (variancePartition).

- The network file must be provided explicitly via the '--decoupler_network' parameter. This file should be in long format and contain at least the source and target columns, with optional weight and sign columns describing the strength and direction of each interaction.

#### Parameters

The Decoupler module includes a min_n parameter to fine-tune its behavior.

- `--decoupler_min_n`: This parameter controls the minimum number of targets a regulator (source) must have in the network to be included in the analysis. Any regulator with fewer than min_n targets will be removed from the network before activity inference is performed.

By default, `--decoupler_min_n` is set to 5, meaning all sources with at least 5 target genes will be evaluated. You can increase this value to filter out poorly supported regulators and reduce noise.

- `--decoupler_methods`: This parameter lets you specify which statistical methods decoupler will use to estimate regulator activities. Decoupler supports multiple methods, each using a different algorithm or statistical approach. You can specify one or more methods by passing them as a comma-separated list.

Example: `--decoupler_methods` mlm,ulsm

##### Network Sources

You can obtain regulatory networks from well-established databases and tools. Common examples include:

- DoRothEA – transcription factor-target interactions (TFs) [DoRothEA](https://www.bioconductor.org/packages/release/data/experiment/html/dorothea.html)

- CollecTRI – curated transcriptional regulatory interactions (TFs) [CollectTRI](https://github.com/saezlab/CollecTRI)

- PROGENy – pathway-responsive gene signatures (pathways) [PROGENy](https://saezlab.github.io/progeny/)

If you want to see the full list of available methods and functions, refer to the function's [official guide](https://decoupler-py.readthedocs.io/en/latest/notebooks/benchmark.html#Multiple-networks).

**Note**: These resources mentioned above are provided only for human or mouse datasets. Please ensure your organism is compatible before enabling this module or provide a custom, species-specific dataset.

### grea

[grea](https://github.com/tpq/propr) (Graph Reverse-Engineered Enrichment Analysis) is the gene set enrichment counterpart to the `propd` differential method, implemented in the [propr](https://github.com/tpq/propr) R package. Unlike GSEA / gProfiler2 / Decoupler, which consume per-gene differential statistics, grea consumes the gene-by-gene **adjacency matrix** produced by propd, so it can only be combined with `--differential_method propd`.

```bash
--differential_method propd \
--functional_method grea \
--gene_sets_files gene_sets.gmt
```

For convenience the pipeline also provides bundled profiles:

```bash
nextflow run nf-core/differentialabundance -profile rnaseq_propd_grea,docker ...
```

Key parameters: `--grea_set_min` and `--grea_set_max` bound gene set sizes; `--grea_permutation` controls the number of permutations used to estimate significance. See the [parameter reference](https://nf-co.re/differentialabundance/parameters) for the full list of `propd_*` and `grea_*` options.

## propd: differential proportionality

[propd](https://github.com/tpq/propr) tests for differential proportionality between groups in compositional data. It is normalization-free (it works on the raw counts via log-ratios) and emits a different set of result columns from DESeq2/Limma/DREAM, most importantly:

- `LFC` - log2 fold-change between groups, used as the differential effect-size column.
- `rcDdis` - reverse cumulative degree distribution; the propd-equivalent of an adjusted p-value, bounded `[0, 1]`. Used as both the p-value and q-value column in pipeline reports.

Because propd does not produce a normalised matrix, downstream plots and reports operate on the raw assay only. The pipeline exposes a `rnaseq_propd` profile for propd-only runs and `rnaseq_propd_grea` to pair propd with grea functional enrichment.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/differentialabundance \
    -profile rnaseq,docker \
    --input samplesheet.csv \
    --contrasts contrasts.yaml \
    --matrix assay_matrix.tsv \
    --gtf mouse.gtf \
    --outdir <OUTDIR> \
    --report_contributors $'Jane Doe\nDirector of Institute of Microbiology\nUniversity of Smallville;John Smith\nPhD student\nInstitute of Microbiology\nUniversity of Smallville'
```

This will launch the pipeline with the `rnaseq` analysis profile (DESeq2 by default) and `docker` for container execution. See below for more information about profiles.

For other data types, use the appropriate analysis profile:

```bash
# Affymetrix microarrays
nextflow run nf-core/differentialabundance \
    -profile affy,docker \
    --input samplesheet.csv \
    --contrasts contrasts.yaml \
    --affy_cel_files_archive cel_files.tar \
    --outdir <OUTDIR>

# MaxQuant proteomics
nextflow run nf-core/differentialabundance \
    -profile maxquant,docker \
    --input samplesheet.csv \
    --contrasts contrasts.yaml \
    --matrix proteinGroups.txt \
    --outdir <OUTDIR>

# GEO SOFT files
nextflow run nf-core/differentialabundance \
    -profile soft,docker \
    --input samplesheet.csv \
    --contrasts contrasts.yaml \
    --querygse GSE12345 \
    --outdir <OUTDIR>

# Generic pre-scaled matrix (Limma)
nextflow run nf-core/differentialabundance \
    -profile generic_matrix,docker \
    --input samplesheet.csv \
    --contrasts contrasts.yaml \
    --matrix my_matrix.tsv \
    --outdir <OUTDIR>

```

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Hints and tips

- If you don't like the colors used in the report, try a different `RColorBrewer` palette by changing the `exploratory_palette_name` and/or `differential_palette_name` parameters.
- In rare cases, some users have reported issues with DESeq2 using all available cores on a machine, rather than those specified in the process configuration. This can be prevented by setting the `OPENBLAS_NUM_THREADS` environment variable.
- By default, `--round_digits` is disabled (`-1`) to avoid unintentional information loss in small numeric values. Enable it only when you explicitly want rounded report tables.

### Scaling up to large sample numbers

#### Deactivating reporting processes

A number of workflow steps are not optimised to deal with large sample numbers and will cause the overall workflow to fail. If you have sample numbers on the order of 100s or more, you should disable these processes like:

```
process {
    withName:'PLOT_EXPLORATORY|PLOT_DIFFERENTIAL|QUARTONOTEBOOK|MAKE_REPORT_BUNDLE|SHINYNGS_APP'{
        ext.when = false
    }
}
```

You will not get the final reporting outcomes of the workflow, but you will get the differential tables produced by DESeq2, Limma, or DREAM, and the results of any gene sets analysis you have enabled.

We have also added a dedicated pipeline parameter, `--skip_reports` that allows you to skip only the Quarto notebook and bundled report while leaving other reporting processes active. The `QUARTONOTEBOOK` process assumes that every grouping variable you pass to it (from the contrasts file’s variable column or PCA-derived informative_variables) exists as a valid, named column in your sample metadata. If you know your metadata or contrasts might be incomplete or non-standard (such as using formula-based yaml files), the you can use this flag to skip these steps.

#### Restricting samples considered by DESeq2, Limma, or DREAM

By default, the DESeq2, Limma, or DREAM differential modules model all samples at once, rather than just the samples involved in the contrast. This is usually the correct thing to do, but when there are large numbers of samples involved in each contrast it may be unnecessary, and things can be sped up significantly by setting `--differential_subset_to_contrast_samples`. This will remove any samples not relevant to the contrast before the main differential analysis routines are called.

### Params files

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/differentialabundance -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
report_contributors: |
  Jane Doe
  Director of Institute of Microbiology
  University of Smallville;John Smith
  PhD student
  Institute of Microbiology
  University of Smallville
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/differentialabundance
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/differentialabundance releases page](https://github.com/nf-core/differentialabundance/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducibility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile rnaseq,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

#### Analysis profiles

These profiles configure the pipeline for specific study types, differential methods, and optional functional enrichment. Combine one analysis profile with a container profile (e.g. `-profile rnaseq,docker`).

**RNA-seq profiles:**

- `rnaseq` — RNA-seq with DESeq2 (default)
- `rnaseq_deseq2_gsea` — RNA-seq with DESeq2 + GSEA enrichment
- `rnaseq_deseq2_gprofiler2` — RNA-seq with DESeq2 + g:Profiler2 enrichment
- `rnaseq_limma` — RNA-seq with Limma (voom)
- `rnaseq_limma_gsea` — RNA-seq with Limma + GSEA enrichment
- `rnaseq_limma_gprofiler2` — RNA-seq with Limma + g:Profiler2 enrichment
- `rnaseq_limma_decoupler` — RNA-seq with Limma + Decoupler enrichment
- `rnaseq_dream` — RNA-seq with DREAM (mixed-effects models)
- `rnaseq_dream_decoupler` — RNA-seq with DREAM + Decoupler enrichment

**Affymetrix profiles:**

- `affy` — Affymetrix microarrays with Limma
- `affy_limma_gsea` — Affymetrix with Limma + GSEA enrichment
- `affy_limma_gprofiler2` — Affymetrix with Limma + g:Profiler2 enrichment

**Other profiles:**

- `maxquant` — MaxQuant proteomics with Limma
- `soft` — GEO SOFT files with Limma
- `generic_matrix` — Generic pre-scaled matrix with Limma
- `generic_matrix_dream` — Generic pre-scaled matrix with DREAM (mixed-effects models)

> [!WARNING]
> Do not override `--differential_method` on the command line when using an analysis profile. Each profile also sets method-specific parameters (e.g. logFC and p-value column names). To change methods, use the appropriate profile (e.g. `-profile rnaseq_limma`).

#### Container profiles

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

#### Other profiles

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `cache_default`
  - The pipeline uses by default `cache = deep` for certain processes downstream `VALIDATOR`, which allows more efficient resuming by checking for input file content to be included in the cache keys. The profile `cache_default` provides users with a the default cache profile (`cache = true`) used in most Nextflow pipelines, in case the `deep` configuration causes issues. See the [Nextflow cache documentation](https://www.nextflow.io/docs/latest/reference/process.html#process-cache) for more information.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/differentialabundance/blob/dev/conf/base.config#L17) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/differentialabundance pipeline is failing after multiple re-submissions of the `DESEQ2_DIFFERENTIAL` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_DIFFERENTIALABUNDANCE:DIFFERENTIALABUNDANCE::DESEQ2_DIFFERENTIAL ([variable:treatment, reference:WT, target:P23H, blocking:, id:treatment_WT_P23H_)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_DIFFERENTIALABUNDANCE:DIFFERENTIALABUNDANCE::DESEQ2_DIFFERENTIAL ([variable:treatment, reference:WT, target:P23H, blocking:, id:treatment_WT_P23H_)'

Caused by:
    Process `NFCORE_DIFFERENTIALABUNDANCE:DIFFERENTIALABUNDANCE::DESEQ2_DIFFERENTIAL (WT_REP1)` terminated with an error exit status (137)

Command executed:
    template 'deseq_de.R'

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    template 'deseq_de.R'
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
