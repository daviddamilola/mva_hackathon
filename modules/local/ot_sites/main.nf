process OT_SITES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    tuple val(meta), path(stage1a_vcf), path(stage1a_csi)
    path ot_jsons

    output:
    tuple val(meta), path("${meta.id}.ot_sites.tsv"), path("${meta.id}.ot_regions.bed"), emit: sites
    path "${meta.id}.ot_sites.summary.json", emit: summary
    path "versions.yml"                    , emit: versions

    script:
    """
    ot_sites.py --genotyped ${stage1a_vcf} --out-tsv ${meta.id}.ot_sites.tsv --out-bed ${meta.id}.ot_regions.bed \\
        ${ot_jsons} > ${meta.id}.ot_sites.summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.ot_sites.tsv ${meta.id}.ot_regions.bed
    echo '{"new_positions": 0}' > ${meta.id}.ot_sites.summary.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.11
    END_VERSIONS
    """
}
