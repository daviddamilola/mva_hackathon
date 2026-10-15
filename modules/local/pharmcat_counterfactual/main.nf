/*
 * Copy-number counterfactual: re-run PharmCAT on the same preprocessed VCF with every duplicated
 * outside call collapsed to one copy (e.g. CYP2D6 *4/*35x2 -> *4/*35), and diff the CPIC/DPWG
 * recommendations. Shows what copy-number detection from the re-aligned reads changed. Skipped
 * (with an explicit 'not applicable' record) when no outside call carries a duplication.
 */
process PHARMCAT_COUNTERFACTUAL {
    tag "$meta.id"
    label 'process_low'

    conda "${projectDir}/modules/local/pharmcat_run/environment.yml"   // same PharmCAT build and data as PHARMCAT_RUN
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/pharmcat3:3.4.0--py313h106432d_0' :
        'quay.io/biocontainers/pharmcat3:3.4.0--py313h106432d_0' }"

    input:
    tuple val(meta), path(report_json), path(preprocessed), path(outside_calls, stageAs: 'outside_??.tsv')

    output:
    path "${meta.id}.pgx_counterfactual.json", emit: json
    path "versions.yml"                      , emit: versions

    script:
    def script_md5 = file("${projectDir}/bin/pgx_counterfactual.py").text.md5()
    """
    cat outside_*.tsv | awk 'NF' > outside_calls.tsv
    pgx_counterfactual.py make outside_calls.tsv counterfactual_calls.tsv changes.json
    if [ "\$(cat changes.json)" != "{}" ]; then
        pharmcat -vcf ${preprocessed} -po counterfactual_calls.tsv -o cf -bf ${meta.id} \\
            -reporterJson -rs ${params.pharmcat_reporter_sources}
        pgx_counterfactual.py diff ${report_json} cf/${meta.id}.report.json changes.json ${meta.id}.pgx_counterfactual.json
    else
        pgx_counterfactual.py diff ${report_json} ${report_json} changes.json ${meta.id}.pgx_counterfactual.json
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: \$(pharmcat -version 2>&1 | head -1 | sed 's/PharmCAT //')
        pgx_counterfactual_py_md5: "${script_md5}"
    END_VERSIONS
    """

    stub:
    """
    echo '{"applicable": false, "reason": "stub"}' > ${meta.id}.pgx_counterfactual.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: 3.4.0
    END_VERSIONS
    """
}
