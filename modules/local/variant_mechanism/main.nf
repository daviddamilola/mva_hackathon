/*
 * Variant-level mechanism for the causal gene, derived rather than asserted:
 *   nmd_ptc.py       is each premature termination codon more than the 50-55 nt rule allows from the
 *                    last exon-exon junction of the MANE Select transcript? (NMD sensitivity)
 *   domain_figure.py where do this child's variants fall in the protein's domain organisation, and
 *                    which domains does a truncating allele keep or lose?
 * Domain coordinates come from a literature asset carrying its own citation; variant positions come
 * from the pipeline's own panel output.
 */
process VARIANT_MECHANISM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(vep_vcf), path(rare_functional)
    path mane_annotation
    path domains

    output:
    path "${meta.id}.nmd_ptc.json"      , emit: nmd_json
    path "${meta.id}.nmd_ptc.tsv"       , emit: nmd_tsv
    path "${meta.id}.domain_map.json"   , emit: domain_json
    path "${meta.id}.domain_map.png"    , emit: figure
    path "versions.yml"                 , emit: versions

    script:
    // digests of the bin scripts: Nextflow does not hash bin/, so they go in the task hash and provenance
    def nmd_md5 = file("${projectDir}/bin/nmd_ptc.py").text.md5()
    def dom_md5 = file("${projectDir}/bin/domain_figure.py").text.md5()
    """
    nmd_ptc.py --vcf ${vep_vcf} --annotation ${mane_annotation} \\
        --out-tsv ${meta.id}.nmd_ptc.tsv --out-json ${meta.id}.nmd_ptc.json

    domain_figure.py --domains ${domains} --rare-functional ${rare_functional} \\
        --gene ${params.mechanism_gene} \\
        --out-png ${meta.id}.domain_map.png --out-json ${meta.id}.domain_map.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)')
        domain_source: "\$(grep -m1 '^# Source' ${domains} | sed 's/^# Source: //')"
        nmd_ptc_py_md5: "${nmd_md5}"
        domain_figure_py_md5: "${dom_md5}"
    END_VERSIONS
    """

    stub:
    """
    echo '{"ptcs": []}' > ${meta.id}.nmd_ptc.json
    echo '{"variants": []}' > ${meta.id}.domain_map.json
    touch ${meta.id}.nmd_ptc.tsv ${meta.id}.domain_map.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.11
    END_VERSIONS
    """
}
