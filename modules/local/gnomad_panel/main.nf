/*
 * gnomAD v4.1 joint (exomes + genomes) frequencies for the panel regions only, fetched remotely
 * from the public bucket with tabix-indexed range requests. Only public gene coordinates are
 * sent; no subject data leaves the machine. Output uses no-chr contigs for VEP --custom.
 */
process GNOMAD_PANEL {
    tag "gnomAD v${params.gnomad_version}"
    label 'process_single'
    maxRetries 3
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.22--h3a4d415_1' :
        'quay.io/biocontainers/bcftools:1.22--h3a4d415_1' }"

    input:
    path panel_bed

    output:
    tuple path("gnomad_joint_panel.vcf.gz"), path("gnomad_joint_panel.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def fields = params.gnomad_fields.tokenize(',').collect { 'INFO/' + it.trim() }.join(',')
    """
    set -o pipefail
    grep -v '^#' ${panel_bed} | cut -f1-3 | sort -k1,1V -k2,2n > regions.nochr.bed
    for c in \$(seq 1 22) X Y; do echo "chr\$c \$c"; done > chr_to_nochr.txt
    parts=()
    for c in \$(cut -f1 regions.nochr.bed | sort -u -V); do
        awk -v c=\$c '\$1==c {print "chr"\$0}' regions.nochr.bed > r_\$c.bed
        url="${params.gnomad_base_url}/gnomad.joint.v${params.gnomad_version}.sites.chr\$c.vcf.bgz"
        bcftools view -R r_\$c.bed "\$url" -Ou \\
            | bcftools annotate -x '^${fields}' -Ou \\
            | bcftools annotate --rename-chrs chr_to_nochr.txt -Oz -o part_\$c.vcf.gz
        parts+=(part_\$c.vcf.gz)
    done
    bcftools concat "\${parts[@]}" -Ou | bcftools sort -Oz -o gnomad_joint_panel.vcf.gz
    tabix -p vcf gnomad_joint_panel.vcf.gz
    rm -f part_*.vcf.gz *.tbi.tmp gnomad.joint.*.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | cut -d' ' -f2)
        gnomad_version: "${params.gnomad_version} (joint exomes+genomes)"
        gnomad_source: "${params.gnomad_base_url}"
        gnomad_records_in_panel: "\$(bcftools view -H gnomad_joint_panel.vcf.gz | wc -l | tr -d ' ')"
    END_VERSIONS
    """

    stub:
    """
    touch gnomad_joint_panel.vcf.gz gnomad_joint_panel.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.22
        gnomad_version: "${params.gnomad_version} (joint exomes+genomes)"
        gnomad_source: "${params.gnomad_base_url}"
        gnomad_records_in_panel: "stub"
    END_VERSIONS
    """
}
