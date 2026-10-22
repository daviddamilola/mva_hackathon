/*
 * The drugged neighbourhood of the mitotic machinery: for each gene in a curated set, every clinical
 * candidate against it, that candidate's mechanism action type, and whether the gene carries one of
 * the mitotic-fidelity GO terms the anti-target rule uses. The gene set is a repository asset with
 * its rationale in the file header; the drugs and their directions are fetched here, never asserted.
 * Public reference data only; no subject data is sent to the API.
 */
process OT_NEIGHBOURHOOD {
    label 'process_single'
    maxRetries 3
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    // the same interpreter-only environment as OT_QUERY; pointing at that file rather than a copy
    // keeps both processes on one cached conda prefix instead of building a second identical one
    conda "${projectDir}/modules/local/ot_query/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    path genes_tsv

    output:
    path "ot_neighbourhood.json", emit: json
    path "versions.yml"         , emit: versions

    script:
    """
    ot_query.py --out ot_neighbourhood.json neighbourhood --genes-tsv ${genes_tsv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        open_targets_endpoint: "https://api.platform.opentargets.org/api/v4/graphql"
        mitotic_drug_targets: "\$(grep -vc '^#' ${genes_tsv}) rows incl. header"
    END_VERSIONS
    """

    stub:
    """
    echo '{"genes": [], "summary": {"genes": 0, "genes_with_candidates": 0, "mechanism_rows": 0, "inhibitory_rows": 0, "genes_without_candidates": [], "genes_without_mitotic_go": [], "non_inhibitory": []}}' > ot_neighbourhood.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.11
        open_targets_endpoint: "https://api.platform.opentargets.org/api/v4/graphql"
    END_VERSIONS
    """
}
