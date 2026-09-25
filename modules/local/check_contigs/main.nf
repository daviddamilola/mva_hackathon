/*
 * Fail fast if the rebuilt reference does not match the subject VCF's contig names and lengths.
 * A mismatch would silently shift every coordinate-based analysis. Reads VCF header only.
 */
process CHECK_CONTIGS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    tuple val(meta), path(fai)
    path vcf

    output:
    path "contig_check.json", emit: report
    path "versions.yml"     , emit: versions

    script:
    """
    check_contigs.py ${vcf} ${fai} ${meta.id} contig_check.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    check_contigs.py ${vcf} ${fai} ${meta.id} contig_check.json --no-strict

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
