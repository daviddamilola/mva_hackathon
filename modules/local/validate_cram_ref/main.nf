/*
 * First node of every CRAM-consuming branch: prove the CRAM was encoded against this reference.
 * (1) every @SQ SN/LN/M5 in the CRAM header must equal the reference dict;
 * (2) a small region must decode with the explicit reference.
 * Downstream processes take this task's output as an input, so they cannot start if it fails.
 */
process VALIDATE_CRAM_REF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'quay.io/biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(cram), path(crai)
    tuple val(ref_meta), path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path(cram), path(crai), path("${meta.id}.cram_ref_check.json"), emit: cram
    path "versions.yml"                                                         , emit: versions

    script:
    """
    samtools head ${cram} > cram_header.sam
    validate_cram_m5.py cram_header.sam ${dict} ${meta.id}.cram_ref_check.json

    # decode test with the explicit reference (fails on reference mismatch or missing blocks)
    samtools view --reference ${fasta} -c ${cram} 1:1000000-1010000 > /dev/null

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    echo '{"pass": true, "stub": true}' > ${meta.id}.cram_ref_check.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.22.1
        python: 3.12.11
    END_VERSIONS
    """
}
