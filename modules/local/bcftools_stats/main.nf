process BCFTOOLS_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.22--h3a4d415_1' :
        'quay.io/biocontainers/bcftools:1.22--h3a4d415_1' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}.bcftools.stats"), emit: stats
    path "versions.yml", emit: versions

    script:
    """
    bcftools stats -f PASS ${vcf} > ${meta.id}.bcftools.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.bcftools.stats
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.22
    END_VERSIONS
    """
}
