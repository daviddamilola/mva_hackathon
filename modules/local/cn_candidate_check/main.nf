/*
 * Read-level check of each exon-run CN candidate from GENE_CN: depth step across the run
 * (including the intron between its exons) vs flanking sequence, and a matching Delly call.
 * A candidate with no step and no Delly call is resolved as noise; anything else stays open.
 */
process CN_CANDIDATE_CHECK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(cram), path(crai), path(ref_check), path(gene_cn)
    path exon_cn
    path sv_panel
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    path "${meta.id}.panel.cn_candidates.tsv" , emit: tsv
    path "${meta.id}.panel.cn_candidates.json", emit: json
    path "versions.yml"                       , emit: versions

    script:
    // digest of the bin script: part of the task hash (Nextflow does not hash bin/), and recorded in provenance
    def script_md5 = file("${projectDir}/bin/cn_candidate_check.py").text.md5()
    """
    cn_candidate_check.py --gene-cn ${gene_cn} --exon-cn ${exon_cn} --sv-panel ${sv_panel} \\
        --cram ${cram} --fasta ${fasta} \\
        --out-tsv ${meta.id}.panel.cn_candidates.tsv --out-json ${meta.id}.panel.cn_candidates.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
        python: \$(python3 --version | cut -d' ' -f2)
        cn_candidate_check_py_md5: "${script_md5}"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.panel.cn_candidates.tsv
    echo '{"stub": true}' > ${meta.id}.panel.cn_candidates.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.22.1
    END_VERSIONS
    """
}
