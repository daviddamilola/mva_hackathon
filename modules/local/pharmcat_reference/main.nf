/*
 * chr-named, bgzipped copy of the primary contigs of *our* reference, for the PharmCAT
 * preprocessor's normalisation step (-refFna). Using the pipeline reference, rather than letting
 * PharmCAT download its own, keeps a single explicit reference for every step.
 */
process PHARMCAT_REFERENCE {
    tag "$ref_meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'quay.io/biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    tuple path("reference.chr.fna.bgz"), path("reference.chr.fna.bgz.fai"), path("reference.chr.fna.bgz.gzi"), emit: fasta
    path "versions.yml", emit: versions

    script:
    """
    set -o pipefail
    samtools faidx ${fasta} \$(seq 1 22) X Y M | sed 's/^>/>chr/' | bgzip -@ ${task.cpus} -c > reference.chr.fna.bgz
    samtools faidx reference.chr.fna.bgz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
        source_reference: "${ref_meta.id}"
    END_VERSIONS
    """

    stub:
    """
    touch reference.chr.fna.bgz reference.chr.fna.bgz.fai reference.chr.fna.bgz.gzi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.22.1
        source_reference: "${ref_meta.id}"
    END_VERSIONS
    """
}
