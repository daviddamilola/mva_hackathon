process QC_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    tuple val(meta), path(samtools_stats), path(flagstat), path(somalier_anc), path(somalier_samples), path(selfsm), path(roh)
    path markdup
    path mosdepth_summary
    path aneuploidy_summary

    output:
    tuple val(meta), path("genome_qc.json"), emit: json
    path "versions.yml", emit: versions

    script:
    """
    qc_summary.py --flagstat ${flagstat} --markdup ${markdup} --mosdepth ${mosdepth_summary} \\
        --somalier-ancestry ${somalier_anc} --somalier-samples ${somalier_samples} --selfsm ${selfsm} \\
        --roh ${roh} --aneuploidy ${aneuploidy_summary} --sample ${meta.id} --declared-sex ${meta.sex} --out genome_qc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    echo '{"checks": []}' > genome_qc.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.11
    END_VERSIONS
    """
}
