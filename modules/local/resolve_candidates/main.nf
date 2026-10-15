/*
 * Resolve the candidate drug list to ChEMBL IDs through the Open Targets search API and apply the
 * anti-target rules. Names in, IDs out: no ChEMBL ID is ever hardcoded in the candidate file, and a
 * name that does not resolve is reported with a status rather than dropped. Excluded candidates stay
 * in the tables with their exclusion reason. Public reference data only; no subject data is sent.
 */
process RESOLVE_CANDIDATES {
    tag "${meta.id}"
    label 'process_single'
    maxRetries 3
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(drug_list)

    output:
    tuple val(meta), path("${meta.id}.track2.candidates.json"), emit: json
    path "${meta.id}.track2.query_ids.txt"                    , emit: query_ids
    path "${meta.id}.track2.resolved.tsv"                     , emit: resolved
    path "${meta.id}.track2.anti_target.tsv"                  , emit: anti_target
    path "versions.yml"                                       , emit: versions

    script:
    def script_md5 = file("${projectDir}/bin/resolve_candidates.py").text.md5()
    """
    resolve_candidates.py --candidates ${drug_list} --prefix ${meta.id}.track2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        open_targets_endpoint: "https://api.platform.opentargets.org/api/v4/graphql"
        drug_list_sha256: "\$(python3 -c 'import hashlib,sys; print(hashlib.sha256(open("${drug_list}","rb").read()).hexdigest())')"
        resolve_candidates_py_md5: "${script_md5}"
    END_VERSIONS
    """

    stub:
    """
    echo '{"candidates": {}, "counts": {}}' > ${meta.id}.track2.candidates.json
    touch ${meta.id}.track2.query_ids.txt ${meta.id}.track2.resolved.tsv ${meta.id}.track2.anti_target.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.11
    END_VERSIONS
    """
}
