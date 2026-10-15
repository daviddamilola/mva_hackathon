process SAMTOOLS_STATS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'quay.io/biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(cram), path(crai), path(ref_check)
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.samtools.stats"), path("${meta.id}.flagstat.txt"), emit: stats
    path "versions.yml", emit: versions

    script:
    """
    samtools stats -@ ${task.cpus} --reference ${fasta} ${cram} > ${meta.id}.samtools.stats
    samtools flagstat -@ ${task.cpus} --input-fmt-option reference=${fasta} ${cram} > ${meta.id}.flagstat.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.samtools.stats ${meta.id}.flagstat.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.22.1
    END_VERSIONS
    """
}
