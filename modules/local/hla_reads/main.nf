/*
 * Reads for HLA typing. The reference has no HLA/ALT contigs, so HLA-derived reads either map to the
 * MHC on chr6 or stay unmapped (divergent alleles). Extracted: the MHC region plus unplaced unmapped
 * reads (index-jumped via '*'), mates re-paired with collate, singletons dropped. CRAM decoded with
 * the explicit reference.
 */
process HLA_READS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'quay.io/biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(cram), path(crai), path(ref_check)
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.hla_R1.fq.gz"), path("${meta.id}.hla_R2.fq.gz"), emit: fastq
    path "versions.yml", emit: versions

    script:
    """
    set -o pipefail
    samtools view -@ ${task.cpus} --reference ${fasta} -u ${cram} ${params.hla_region} '*' \\
        | samtools collate -@ ${task.cpus} -u -O -T collate_tmp - \\
        | samtools fastq -@ ${task.cpus} -n -1 ${meta.id}.hla_R1.fq.gz -2 ${meta.id}.hla_R2.fq.gz -0 /dev/null -s /dev/null -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
        hla_region: "${params.hla_region} + unplaced unmapped"
        read_pairs: "\$(( \$(gzip -dc ${meta.id}.hla_R1.fq.gz | wc -l) / 4 ))"
    END_VERSIONS
    """

    stub:
    """
    echo | gzip -c > ${meta.id}.hla_R1.fq.gz; echo | gzip -c > ${meta.id}.hla_R2.fq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.22.1
        hla_region: "${params.hla_region} + unplaced unmapped"
    END_VERSIONS
    """
}
