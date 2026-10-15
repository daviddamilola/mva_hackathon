process ALDY_PARSE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    tuple val(meta), path(aldy_files)

    output:
    tuple val(meta), path("${meta.id}.aldy_calls.json"), path("${meta.id}.outside_calls.tsv"), emit: calls
    path "versions.yml", emit: versions

    script:
    """
    aldy_parse.py --outside-genes '${params.pgx_outside_call_genes}' \\
        --json ${meta.id}.aldy_calls.json --outside ${meta.id}.outside_calls.tsv ${aldy_files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        outside_call_genes: "${params.pgx_outside_call_genes}"
    END_VERSIONS
    """

    stub:
    """
    echo '{}' > ${meta.id}.aldy_calls.json
    touch ${meta.id}.outside_calls.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.11
        outside_call_genes: "${params.pgx_outside_call_genes}"
    END_VERSIONS
    """
}
