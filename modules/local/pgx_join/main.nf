process PGX_JOIN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    tuple val(meta), path(calls_json), path(pharmcat_report), path(stage1a_vcf), path(ot_vcf)
    path drug_jsons
    path target_json
    path meta_json
    path candidates_json
    val drug_list_label
    val drug_list_example

    output:
    tuple val(meta), path("${meta.id}.pgx_ot.verdicts.tsv"), emit: verdicts
    path "${meta.id}.pgx_ot.forward_annotations.tsv"       , emit: forward
    path "${meta.id}.pgx_ot.evidence_levels.tsv"           , emit: levels
    path "${meta.id}.pgx_ot.summary.json"                  , emit: summary
    path "versions.yml"                                    , emit: versions

    script:
    def direct = params.ot_include_direct_targets ? '--include-direct-targets' : ''
    def example = drug_list_example ? '--drug-list-example' : ''
    """
    pgx_join.py --calls ${calls_json} --pharmcat-report ${pharmcat_report} \\
        --vcf ${stage1a_vcf} ${ot_vcf} --drugs ${drug_jsons} --targets ${target_json} --meta ${meta_json} \\
        --candidates ${candidates_json} --assert-levels '${params.ot_assert_levels}' ${direct} \\
        --drug-list-label '${drug_list_label}' ${example} --prefix ${meta.id}.pgx_ot

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        open_targets_release: "\$(python3 -c 'import json; print(json.load(open("${meta_json}"))["data_release"])')"
        open_targets_api: "\$(python3 -c 'import json; print(json.load(open("${meta_json}"))["api_version"])')"
        asserted_levels: "${params.ot_assert_levels}"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.pgx_ot.verdicts.tsv ${meta.id}.pgx_ot.forward_annotations.tsv ${meta.id}.pgx_ot.evidence_levels.tsv
    echo '{"stub": true}' > ${meta.id}.pgx_ot.summary.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.11
        open_targets_release: "stub"
        open_targets_api: "stub"
        asserted_levels: "${params.ot_assert_levels}"
    END_VERSIONS
    """
}
