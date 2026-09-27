process MULTIQC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.35--pyhdfd78af_2' :
        'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_2' }"

    input:
    path qc_files, stageAs: 'qc/*'

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data"       , emit: data
    path "versions.yml"       , emit: versions

    script:
    """
    multiqc --force --title 'mva-track2 QC' --filename multiqc_report.html qc/
    mv multiqc_report_data multiqc_data

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/.*version //')
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html; mkdir multiqc_data
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: 1.35
    END_VERSIONS
    """
}
