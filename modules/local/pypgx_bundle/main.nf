/*
 * pypgx-bundle resource directory (1KGP reference haplotype panel, CNV models), pinned to the
 * PyPGx version and verified against a recorded sha256. storeDir: downloaded once, then reused.
 * The bundle repository states no licence; this is recorded in provenance.
 */
process PYPGX_BUNDLE {
    tag "pypgx-bundle-${version}"
    label 'process_single'
    storeDir "${params.reference_dir}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    tuple val(version), val(sha256)

    output:
    path "pypgx-bundle-${version}", emit: bundle

    script:
    """
    set -o pipefail
    curl -fsSL --retry 5 -o bundle.tar.gz https://github.com/sbslee/pypgx-bundle/archive/refs/tags/${version}.tar.gz
    verify_checksum.py sha256 ${sha256} bundle.tar.gz
    tar -xzf bundle.tar.gz
    rm bundle.tar.gz
    """

    stub:
    """
    mkdir -p pypgx-bundle-${version}
    """
}
