# nf-core/differentialabundance: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/differentialabundance/usage](https://nf-co.re/differentialabundance/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

Differential analysis is a common task in a variety of use cases. In essence, all these use cases entail taking an input matrix containing features (e.g. genes) and observations (e.g. samples), and comparing groups of observations in all or a subset of the features. The feature/ observation language here reflects our hope that this workflow will extend in future to encompass a variety of applications where an assumption of gene vs sample may not be a valid one- though that is the application to which the first release will apply.

With the above in mind, running this workflow requires:

- a matrix of quantifications with observations by column and features by row
- a description of the observations such as a sample sheet from RNA-seq analysis
- a description of the features, for our initial RNA-seq application this is simply the GTF file from which gene annotations can be derived
- a specification of how the matrix should be split, and how the resulting groups should be compared

Currently, this limits the usage of differentialabundance to genomics applications where a GTF annotation file is available. Future developments of the pipeline will encompass changes that enable users to use the pipeline on more generic cases.

## Observations (samplesheet) input

```bash
--input '[path to samplesheet file]'
```

This may well be the same sample sheet used to generate the input matrix. For example, in RNA-seq this might be the same sample sheet, perhaps derived from [fetchngs](https://github.com/nf-core/fetchngs), that was input to the [RNA-seq workflow](https://github.com/nf-core/rnaseq). It may be necessary to add columns that describe the groups you want to compare.

For example:

```console
sample,fastq_1,fastq_2,condition,replicate,batch
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,control,1,A
CONTROL_REP2,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,control,2,B
CONTROL_REP3,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,control,3,A
TREATED_REP1,AEG588A2_S1_L002_R1_001.fastq.gz,AEG588A2_S1_L002_R2_001.fastq.gz,treated,1,B
TREATED_REP2,AEG588A2_S1_L003_R1_001.fastq.gz,AEG588A2_S1_L003_R2_001.fastq.gz,treated,2,A
TREATED_REP3,AEG588A2_S1_L004_R1_001.fastq.gz,AEG588A2_S1_L004_R2_001.fastq.gz,treated,3,B
```

The file can be tab or comma separated.

## Matrix file

```bash
--matrix '[path to matrix file]'
```

This is a numeric square matrix file, comma or tab-separated, with a column for every observation, and features corresponding to the supplied feature set. The parameters `--observations_id_col` and `--features_id_col` define which of the associated fields should be matched in those inputs.

## Contrasts file

```bash
--contrasts '[path to contrasts file]'
```

The contrasts file references the observations file to define groups of samples to compare. For example, based on the sample sheet above we could define contrasts like:

```console
id,variable,reference,target,blocking
condition_control_treated,condition,control,treated,
condition_control_treated_blockrep,condition,control,treated,replicate;batch
```

The necessary fields in order are:

- `id` - an arbitrary identifier, will be used to name contrast-wise output files
- `variable` - which column from the observations information will be used to define groups
- `reference` - the base/ reference level for the comparison. If features have higher values in this group than target they will generate negative fold changes
- `target` - the target/ non-reference level for the comparison. If features have higher values in this group than the reference they will generate positive fold changes
- `blocking` - semicolon-delimited, any additional variables (also observation columns) that should be modelled alongside the contrast variable

The file can be tab or comma separated.

## GTF file

```bash
--gtf '[path to gtf file]'
```

This is the format currently used to supply feature metadata. In the RNA-seq case it should match the GTF used in that workflow.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/differentialabundance \
    --input samplesheet.csv \
    --contrasts contrasts.csv \
    --matrix assay_matrix.tsv \
    --gtf mouse.gtf \
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

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/differentialabundance
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/differentialabundance releases page](https://github.com/nf-core/differentialabundance/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

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
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.

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

#### For beginners

A first step to bypass this error, you could try to increase the amount of CPUs, memory, and time for the whole pipeline. Therefor you can try to increase the resource for the parameters `--max_cpus`, `--max_memory`, and `--max_time`. Based on the error above, you have to increase the amount of memory. Therefore you can go to the [parameter documentation of differentialabundance](https://nf-co.re/differentialabundance/1.0.0/parameters) and scroll down to the `show hidden parameter` button to get the default value for `--max_memory`. In this case 128GB, you than can try to run your pipeline again with `--max_memory 200GB -resume` to skip all process, that were already calculated. If you can not increase the resource of the complete pipeline, you can try to adapt the resource for a single process as mentioned below.

#### Advanced option on process level

To bypass this error you would need to find exactly which resources are set by the `DESEQ2_DIFFERENTIAL` process. The quickest way is to search for `process DESEQ2_DIFFERENTIAL` in the [nf-core/differentialabundance Github repo](https://github.com/nf-core/differentialabundance/search?q=process+DESEQ2_DIFFERENTIAL).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/deseq2/differential/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_medium`](https://github.com/nf-core/differentialabundance/blob/0c0e457d88a33baeb3061114978eebee6cceb46c/modules/nf-core/deseq2/differential/main.nf#L3).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_medium` label are set in the pipeline's [`base.config`](https://github.com/nf-core/differentialabundance/blob/0c0e457d88a33baeb3061114978eebee6cceb46c/conf/base.config#L37) which in this case is defined as 2GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `DESEQ2_DIFFERENTIAL` process failure by creating a custom config file that sets a higher memory.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_DIFFERENTIALABUNDANCE:DIFFERENTIALABUNDANCE::DESEQ2_DIFFERENTIAL' {
        memory = 20.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_DIFFERENTIALABUNDANCE:DIFFERENTIALABUNDANCE::DESEQ2_DIFFERENTIAL` in the config file because this takes priority over the short name (`DESEQ2_DIFFERENTIAL`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers (advanced users)

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

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
