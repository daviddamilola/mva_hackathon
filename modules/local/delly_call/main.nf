/*
 * Germline structural variants genome-wide with Delly 2.x 'sr' (short-read; paired-end + split-read
 * evidence; DEL, DUP, INV, INS, BND), CRAM decoded with the explicit reference, Delly's hg38
 * exclude list (telomeres/centromeres/low-mappability) converted to no-chr names. Genome-wide rather
 * than panel-restricted, so translocations with a partner outside the panel are still seen.
 */
process DELLY_CALL {
    tag "$meta.id"
    label 'process_medium'
    label 'process_long'
    memory { 8.GB * task.attempt }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/delly:2.6.0--h3752d28_0' :
        'quay.io/biocontainers/delly:2.6.0--h3752d28_0' }"

    input:
    tuple val(meta), path(cram), path(crai), path(ref_check)
    tuple val(ref_meta), path(fasta), path(fai)
    tuple val(excl_name), path(exclude)

    output:
    tuple val(meta), path("${meta.id}.delly.vcf.gz"), path("${meta.id}.delly.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    delly sr -g ${fasta} -x ${exclude} -q ${params.delly_min_mapq} -h ${task.cpus} -o ${meta.id}.delly.bcf ${cram}
    bcftools view ${meta.id}.delly.bcf -Oz -o ${meta.id}.delly.vcf.gz
    tabix -p vcf ${meta.id}.delly.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
        bcftools: \$(bcftools --version | head -1 | cut -d' ' -f2)
        delly_args: "sr -g <reference> -x <hg38 exclude, no-chr> -q ${params.delly_min_mapq}"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.delly.vcf.gz ${meta.id}.delly.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: 2.6.0
        bcftools: 1.24
        delly_args: "sr -g <reference> -x <hg38 exclude, no-chr> -q ${params.delly_min_mapq}"
    END_VERSIONS
    """
}
