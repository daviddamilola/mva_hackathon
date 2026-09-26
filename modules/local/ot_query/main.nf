/*
 * Open Targets Platform GraphQL queries (public reference data only; no subject data is sent).
 * mode: 'meta' | 'drug' (id = ChEMBL ID) | 'target' (id = comma-separated gene symbols)
 */
process OT_QUERY {
    tag "${mode}:${id}"
    label 'process_single'
    maxRetries 3
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    tuple val(mode), val(id)

    output:
    tuple val(mode), val(id), path("ot_${mode}_${safe}.json"), emit: json
    path "versions.yml"                                     , emit: versions

    script:
    safe = id.replaceAll(/[^A-Za-z0-9_-]/, '_').take(60)
    def args = mode == 'drug' ? "drug --chembl ${id}" : (mode == 'target' ? "target --genes '${id}'" : 'meta')
    """
    ot_query.py --out ot_${mode}_${safe}.json ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        open_targets_endpoint: "https://api.platform.opentargets.org/api/v4/graphql"
    END_VERSIONS
    """

    stub:
    safe = id.replaceAll(/[^A-Za-z0-9_-]/, '_').take(60)
    """
    if [ "${mode}" = "meta" ]; then
        echo '{"api_version": "stub", "data_release": "stub"}' > ot_${mode}_${safe}.json
    elif [ "${mode}" = "drug" ]; then
        echo '{"query_id": "${id}", "found": true, "resolution": [], "rows": []}' > ot_${mode}_${safe}.json
    else
        echo '{"genes": {}}' > ot_${mode}_${safe}.json
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.11
        open_targets_endpoint: "https://api.platform.opentargets.org/api/v4/graphql"
    END_VERSIONS
    """
}
