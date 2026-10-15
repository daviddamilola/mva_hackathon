/*
 * Aldy vs PyPGx (Stargazer successor) star-allele agreement, per gene. A report, not a gate:
 * differences are listed, not failed on.
 */
process STAR_ALLELE_COMPARE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    tuple val(meta), path(aldy_json), path(pypgx_tsvs)

    output:
    tuple val(meta), path("${meta.id}.aldy_vs_pypgx.tsv"), emit: tsv
    path "${meta.id}.pypgx_calls.tsv"                     , emit: pypgx
    path "${meta.id}.aldy_vs_pypgx.summary.json"          , emit: summary
    path "versions.yml"                                   , emit: versions

    script:
    """
    head -1 \$(ls ${pypgx_tsvs} | head -1) > ${meta.id}.pypgx_calls.tsv
    for f in ${pypgx_tsvs}; do tail -n +2 \$f; done | sort >> ${meta.id}.pypgx_calls.tsv
    compare_star_alleles.py --aldy ${aldy_json} --other ${meta.id}.pypgx_calls.tsv \\
        --out ${meta.id}.aldy_vs_pypgx.tsv 2> ${meta.id}.aldy_vs_pypgx.summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.aldy_vs_pypgx.tsv ${meta.id}.pypgx_calls.tsv
    echo '{"stub": true}' > ${meta.id}.aldy_vs_pypgx.summary.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.11
    END_VERSIONS
    """
}
