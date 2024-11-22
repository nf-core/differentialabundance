# nf-core/differentialabundance: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/differentialabundance/usage](https://nf-co.re/differentialabundance/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

Differential analysis is a common task in a variety of use cases. In essence, all these use cases entail taking an input matrix containing features (e.g. genes) and observations (e.g. samples), and comparing groups of observations in all or a subset of the features. The feature/ observation language here reflects our hope that this workflow will extend in future to encompass a variety of applications where an assumption of gene vs sample may not be a valid one- though that is the application to which the first release will apply.

With the above in mind, running this workflow requires:

- a set of abundance values. This can be:
  - (for RNA-seq or MaxQuant proteomics measurements): a matrix of quantifications with observations by column and features by row
  - (for Affymetrix microarrays): a tar'd archive of CEL files
- a description of the observations such as a sample sheet from RNA-seq analysis
- a description of the features, for our initial RNA-seq application this can be simply the GTF file from which gene annotations can be derived. For Affymetrix arrays this can be derived from the array platform annotation package automatically. Skip for MaxQuant. You can also supply your own table.
- a specification of how the matrix should be split, and how the resulting groups should be compared

## Observations (samplesheet) input

```bash
--input '[path to samplesheet file]'
```

This may well be the same sample sheet used to generate the input matrix. For example, in RNA-seq this might be the same sample sheet, perhaps derived from [fetchngs](https://github.com/nf-core/fetchngs), that was input to the [RNA-seq workflow](https://github.com/nf-core/rnaseq). It may be necessary to add columns that describe the groups you want to compare. The columns that the pipeline requires are:

- a column listing the sample IDs (must be the same IDs as in the abundance matrix), in the example below it is called 'sample'. For some study_types, this column might need to be filled in with file names, e.g. when doing an affymetrix analysis.
- one or more columns describing conditions for the differential analysis. In the example below it is called 'condition'
- optionally one or more columns describing sample batches or similar which you want to be considered in the analysis. In the example below it is called 'batch'

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

The file can be tab or comma separated.

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

The "file" column in this example is used to specify the data file associated with each sample, which is essential for data analysis and interpretation.

## Abundance values

### RNA-seq and similar

```bash
--matrix '[path to matrix file]'
```

This is a numeric square matrix file, comma or tab-separated, with a column for every observation, and features corresponding to the supplied feature set. The parameters `--observations_id_col` and `--features_id_col` define which of the associated fields should be matched in those inputs.

#### Outputs from nf-core/rnaseq and other tximport-processed results

The nf-core RNAseq workflow incorporates [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html) for producing quantification matrices. From [version 3.12.2](https://github.com/nf-core/rnaseq/releases/tag/3.13.2), it additionally provides transcript length matrices which can be directly consumed by DESeq2 to model length bias across samples.

To use this approach, include the transcript lengths file with the **raw counts**:

```bash
--matrix 'salmon.merged.gene_counts.tsv' \
--transcript_length_matrix 'salmon.merged.gene_lengths.tsv'
```

Without the transcript lengths, for instance in earlier rnaseq workflow versions, follow the second recommendation in the [tximport documentation](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#Downstream_DGE_in_Bioconductor):

> "Use the tximport argument `countsFromAbundance='lengthScaledTPM'` or `'scaledTPM'`, then employ the gene-level count matrix `txi$counts` directly in downstream software, a method we call 'bias corrected counts without an offset'"

This aligns with the **gene_counts_length_scaled.tsv** or **gene_counts_scaled.tsv** matrices in the rnaseq workflow.

It is important to note that the documentation advises:

> "Do not manually pass the original gene-level counts to downstream methods without an offset."

So we **do not recommend** raw counts files such as `salmon.merged.gene_counts.tsv` as input for this workflow **except** where the transcript lengths are also provided.

### MaxQuant intensities

```bash
--matrix '[path to matrix file]'
```

This is the proteinGroups.txt file produced by MaxQuant. It is a tab-separated matrix file with a column for every observation (plus additional columns for other types of measurements and information); each row contains these data for a set of proteins. The parameters `--observations_id_col` and `--features_id_col` define which of the associated fields should be matched in those inputs. The parameter `--proteus_measurecol_prefix` defines which prefix is used to extract those matrix columns which contain the measurements to be used. For example, the default `LFQ intensity ` will indicate that columns like LFQ intensity S1, LFQ intensity S2, LFQ intensity S3 etc. are used (one whitespace is automatically added if necessary).

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

```bash
--contrasts '[path to contrasts file]'
```

The contrasts file references the observations file to define groups of samples to compare. For example, based on the sample sheet above we could define contrasts like:

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
- `exclude_samples_col` and `exclude_samples_values` - the former being a valid column in the samples sheet, the latter a semicolon-delimited list of values in that column which should be used to select samples prior to differential modelling. This is helpful where certain samples need to be exluded prior to analysis of a given contrast.

The file can be tab or comma separated.

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

## Working with the output R markdown file

The pipeline produces an R markdown file which, if you're proficient in R, you can use to tweak the report after it's generated (**note**- if you need the same customisations repeatedly we would recommend you supply your own template using the `report_file` parameter).

To work with R markdown files you will need Rstudio. You will also need to have the ShinyNGS R module [installed](https://github.com/pinin4fjords/shinyngs#installation), since it supplies a lot of the accessory plotting functions etc that you will need. The exact way you will do this may depend on your exact systems, but for example

### 1. Create a conda environment with Shinyngs and activate it

```bash
conda create -n shinyngs r-shinyngs
conda activate shinyngs
```

### 2. Open RStudio from this environment

For example, on a Mac Terminal:

```bash
open -na Rstudio
```

Now, unzip the report archive, and in RStudio change directory to that location:

```
setwd("/path/to/unzipped/directory")
```

Now open the R Markdown file from the RStudio UI, and you should have everything you need to run the various code segments and render the whole document to HTML again if you wish.

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

Currently, two tools can be used to do gene set enrichment analysis.

### GSEA

[GSEA](https://www.gsea-msigdb.org/gsea/index.jsp) tests for differential genes from within a user-provided set of genes; this requires a GMT or GMX file. The following example shows how to enable this:

```bash
--gsea_run true \
--gene_sets_files gene_sets.gmt
```

### gProfiler2

The [gprofiler2](https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html) package can be used to test which pathways are enriched in the sets of differential genes produced by the the DESeq2 or limma modules. It is an R interface for the gprofiler webtool. In the simplest form, this feature can be enabled with the parameters from the following example:

```bash
--gprofiler2_run true \
--gprofiler2_organism mmusculus
```

If gene sets have been specified to the workflow via `--gene_sets_files` these are used by default. Specifying `--gprofiler2_organism` (mmusculus for Mus musculus, hsapiens for Homo sapiens etc.) will override those gene sets with gprofiler's own for the relevant species. `--gprofiler2_token` will override both options and use gene sets from a previous gprofiler run.

By default the analysis will be run with a background list of genes that passed the abundance filter (i.e. those genes that actually had some expression); see for example https://doi.org/10.1186/s13059-015-0761-7 for why this is advisable. You can provide your own background list with `--gprofiler2_background_file background.txt`or if you want to not use any background, set `--gprofiler2_background_file false`.

Check the [pipeline webpage](https://nf-co.re/differentialabundance/parameters#gprofiler2) for a full listing of the relevant parameters.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/differentialabundance \
    [-profile rnaseq OR -profile rnaseq_limma OR -profile affy] \
    --input samplesheet.csv \
    --contrasts contrasts.csv \
    [--matrix assay_matrix.tsv OR --affy_cel_files_archive cel_files.tar] \
    [--gtf mouse.gtf OR --features features.tsv] \
    --outdir <OUTDIR>  \
    -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

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

### Scaling up to large sample numbers

#### Deactivating reporting processes

A number of workflow steps are not optimised to deal with large sample numbers and will cause the overall workflow to fail. If you have sample numbers on the order of 100s or more, you should disable these processes like:

```
process {
    withName:'PLOT_EXPLORATORY|PLOT_DIFFERENTIAL|RMARKDOWNNOTEBOOK|MAKE_REPORT_BUNDLE|SHINYNGS_APP'{
        ext.when = false
    }
}
```

You will not get the final reporting outcomes of the workflow, but you will get the differential tables produced by DESeq2 or Limma, and the results of any gene seta analysis you have enabled.

#### Restricting samples considered by DESeq2 or Limma

By default, the DESeq2 or Limma differential modules model all samples at once, rather than just the samples involved in the contrast. This is usually the correct thing to do, but when there are are large numbers of samples involved in each contrast it may be unnecessary, and things can be sped up significantly by setting `--differential_subset_to_contrast_samples`. This will remove any samples not relevant to the contrast before the main differential analysis routines are called.

### Params files

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/differentialabundance -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/differentialabundance
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/differentialabundance releases page](https://github.com/nf-core/differentialabundance/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducibility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

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

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

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
