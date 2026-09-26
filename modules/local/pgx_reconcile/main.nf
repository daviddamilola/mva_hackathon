process PGX_RECONCILE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    tuple val(meta), path(phenotype_json), path(aldy_json), path(off_definition)

    output:
    tuple val(meta), path("${meta.id}.pgx_calls.tsv"), path("${meta.id}.pgx_calls.json"), emit: calls
    path "versions.yml", emit: versions

    script:
    """
    pgx_reconcile.py ${phenotype_json} ${aldy_json} ${off_definition} '${params.pgx_genes}' ${meta.id}.pgx_calls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    pgx_reconcile.py ${phenotype_json} ${aldy_json} ${off_definition} '${params.pgx_genes}' ${meta.id}.pgx_calls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.11
    END_VERSIONS
    """
}
