/*
 * Subject variants in the mechanistic panel (gene +/-10 kb), from the distributed VCF (same
 * reference, verified by CHECK_CONTIGS). All FILTER values are kept (reported, not dropped).
 * Multiallelics are split and indels left-normalised against the explicit reference.
 */
process PANEL_VARIANTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.22--h3a4d415_1' :
        'quay.io/biocontainers/bcftools:1.22--h3a4d415_1' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path panel_bed
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.panel.vcf.gz"), path("${meta.id}.panel.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    set -o pipefail
    grep -v '^#' ${panel_bed} | cut -f1-3 | sort -k1,1V -k2,2n > regions.bed
    # the distributed .tbi (old format, no count metadata) supports region queries; never re-index
    # a staged input in place (tabix -f could write through the symlink into data/)
    bcftools view -R regions.bed ${vcf} -Ou \\
        | bcftools norm -m -any -f ${fasta} -Ou \\
        | bcftools view -e 'ALT="*"' -Ou \\
        | bcftools sort -Oz -o ${meta.id}.panel.vcf.gz
    tabix -p vcf ${meta.id}.panel.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.panel.vcf.gz ${meta.id}.panel.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.22
    END_VERSIONS
    """
}
