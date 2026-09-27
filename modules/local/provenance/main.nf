process PROVENANCE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    path versions
    path run_meta
    path pharmcat_report
    path ot_summary

    output:
    path "provenance.json", emit: json

    script:
    def ot = ot_summary.name != 'NO_OT_SUMMARY' ? "--ot-summary ${ot_summary}" : ''
    """
    provenance.py --versions ${versions} --meta ${run_meta} --pharmcat-report ${pharmcat_report} ${ot} --out provenance.json
    """

    stub:
    """
    provenance.py --versions ${versions} --meta ${run_meta} --out provenance.json
    """
}
