process CUSTOM_TABULARTOGSEACLS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5a1deaeed2bbe6deb0312419a847e7005d37f9aaa7ff9f447502ff011fc23ac8/data' :
        'community.wave.seqera.io/library/pandas_python:fd8290c2da2fd6ae' }"

    input:
    tuple val(meta), path(samples)

    output:
    tuple val(meta), path("*.cls"), emit: cls
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: []
    def prefix = task.ext.prefix ?: "${meta.id}"
    def separator = args.separator ?: ( samples.getName().endsWith(".tsv") ? '\t' : ',' )
    separator = separator == '\t' ? '\\t' : separator
    def variable = args.variable
    if ( !variable ) error "Supply a variable in the sample sheet from which to derive classes"

    """
    python3 -c "
    import pandas as pd

    # Parameters
    samplesheet_path = '${samples}'
    cls_file = '${prefix}.cls'
    variable = '${variable}'
    separator = '${separator}'

    # Load the samplesheet
    df = pd.read_csv(samplesheet_path, sep=separator)

    # Check if the variable exists uniquely
    if variable not in df.columns:
        raise ValueError(f'Column \"{variable}\" not found in samplesheet.')
    elif len([col for col in df.columns if col == variable]) > 1:
        raise ValueError(f'Multiple columns matched the variable name \"{variable}\". Ensure unique column names.')

    # Extract classes
    classes = df[variable].fillna('empty').astype(str).tolist()
    unique_classes = list(dict.fromkeys(classes)) # Preserve first occurrence order

    # Write the .cls file
    with open(cls_file, 'w') as f:
        f.write(f'{len(classes)} {len(unique_classes)} 1\\n')
        f.write('#' + ' '.join(unique_classes) + '\\n')
        f.write(' '.join(classes) + '\\n')
    "

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o \"([0-9]{1,}\\.)+[0-9]{1,}\")
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
