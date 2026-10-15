/*
 * Merge coordinate-sorted lane CRAMs (same library) and mark duplicates across lanes.
 * Inputs carry fixmate -m tags (MC, ms) from FETCH_ALIGN_LANE.
 */
process MERGE_MARKDUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'quay.io/biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(crams, stageAs: 'lanes/*')
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.cram"), path("${meta.id}.cram.crai"), emit: cram
    path "${meta.id}.markdup.stats.txt"                                  , emit: stats
    path "versions.yml"                                                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    samtools merge -@ ${task.cpus} -u --reference ${fasta} -o - lanes/*.cram \\
        | samtools markdup -@ ${task.cpus} ${args} --reference ${fasta} \\
            -f ${meta.id}.markdup.stats.txt -T markdup_tmp \\
            -O cram --write-index - ${meta.id}.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
        markdup_args: "${args}"
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch ${meta.id}.cram ${meta.id}.cram.crai ${meta.id}.markdup.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.22.1
        markdup_args: "${args}"
    END_VERSIONS
    """
}
