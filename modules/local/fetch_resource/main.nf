/*
 * Small public resource, pinned by URL (commit-specific where the host is GitHub) and sha256.
 * If the file is gzipped it is decompressed; 'chr' prefixes are stripped to match the no-chr
 * reference; track/header lines are dropped. storeDir: fetched once.
 */
process FETCH_RESOURCE {
    tag "$name"
    label 'process_single'
    storeDir "${params.reference_dir}/resources"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    tuple val(name), val(url), val(sha256)

    output:
    tuple val(name), path("${name}.nochr.txt"), emit: file

    script:
    """
    set -o pipefail
    curl -fsSL --retry 5 -o raw.download '${url}'
    verify_checksum.py sha256 ${sha256} raw.download
    ( gzip -dc raw.download 2>/dev/null || cat raw.download ) \\
        | grep -vE '^(track|browser|#)' \\
        | sed -E 's/^chr//' > ${name}.nochr.txt
    """

    stub:
    """
    touch ${name}.nochr.txt
    """
}
