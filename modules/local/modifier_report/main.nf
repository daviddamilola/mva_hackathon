process MODIFIER_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    tuple val(meta), path(vep_vcf)
    path panel_bed

    output:
    tuple val(meta), path("${meta.id}.modifier_panel.rare_functional.tsv"), emit: rare_functional
    path "${meta.id}.modifier_panel.all.tsv"                               , emit: all
    path "${meta.id}.modifier_panel.summary.json"                          , emit: summary
    path "versions.yml"                                                    , emit: versions

    script:
    """
    modifier_panel_report.py --vcf ${vep_vcf} --bed ${panel_bed} --sample ${meta.id} \\
        --prefix ${meta.id}.modifier_panel --max-af ${params.panel_max_af} --spliceai-min ${params.spliceai_min_delta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        panel_max_af: "${params.panel_max_af}"
        spliceai_min_delta: "${params.spliceai_min_delta}"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.modifier_panel.rare_functional.tsv ${meta.id}.modifier_panel.all.tsv
    echo '{"stub": true}' > ${meta.id}.modifier_panel.summary.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.11
        panel_max_af: "${params.panel_max_af}"
        spliceai_min_delta: "${params.spliceai_min_delta}"
    END_VERSIONS
    """
}
