/*
 * AlphaMissense (hg38) restricted to the panel. Resumable download (the file has no index, so
 * region queries are impossible), MD5 verified against Google Cloud Storage, filtered to the panel,
 * and the full file deleted. Scores are CC BY-NC-SA 4.0 (non-commercial).
 */
process ALPHAMISSENSE_PANEL {
    tag "AlphaMissense"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    path panel_bed

    output:
    tuple path("alphamissense_panel.vcf.gz"), path("alphamissense_panel.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    for i in 1 2 3 4 5 6 7 8; do
        curl -fsSL -C - --retry 5 -o am.tsv.gz ${params.alphamissense_url} && break || sleep 30
    done
    grep -v '^#' ${panel_bed} | cut -f1-3 > panel.bed
    alphamissense_panel.py --input am.tsv.gz --bed panel.bed --md5-b64 '${params.alphamissense_md5_b64}' --out alphamissense_panel.vcf
    rm am.tsv.gz
    bgzip alphamissense_panel.vcf && tabix -p vcf alphamissense_panel.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alphamissense_source: "${params.alphamissense_url}"
        alphamissense_md5_b64: "${params.alphamissense_md5_b64}"
        alphamissense_licence: "CC BY-NC-SA 4.0"
    END_VERSIONS
    """

    stub:
    """
    touch alphamissense_panel.vcf.gz alphamissense_panel.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alphamissense_source: "${params.alphamissense_url}"
        alphamissense_md5_b64: "${params.alphamissense_md5_b64}"
        alphamissense_licence: "CC BY-NC-SA 4.0"
    END_VERSIONS
    """
}
