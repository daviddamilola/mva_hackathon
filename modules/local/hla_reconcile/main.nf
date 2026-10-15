process HLA_RECONCILE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    tuple val(meta), path(t1k), path(optitype)

    output:
    tuple val(meta), path("${meta.id}.hla_outside_calls.tsv"), emit: outside
    path "${meta.id}.hla_calls.tsv"                          , emit: calls
    path "${meta.id}.hla_summary.json"                       , emit: summary
    path "versions.yml"                                      , emit: versions

    script:
    """
    hla_reconcile.py --t1k ${t1k} --optitype ${optitype} --pharmcat-genes '${params.hla_pharmcat_genes}' --prefix ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        hla_resolution: "2-field"
        passed_only_if_concordant: "true"
    END_VERSIONS
    """

    stub:
    """
    hla_reconcile.py --t1k ${t1k} --optitype ${optitype} --pharmcat-genes '${params.hla_pharmcat_genes}' --prefix ${meta.id}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.11
    END_VERSIONS
    """
}
