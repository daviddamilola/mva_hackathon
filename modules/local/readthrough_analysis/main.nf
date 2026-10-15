/*
 * Readthrough characterisation of the causal gene's premature termination codon: which stop codon,
 * what 3' context, and what protein a readthrough agent would actually produce (bin/readthrough_analysis.py).
 *
 * The near-cognate products are reachable by one base change from the MUTANT stop codon but not from
 * the reference codon, so AlphaMissense's reference-anchored SNV file (used elsewhere in this pipeline)
 * cannot score them. This streams the all-amino-acid-substitution file and keeps one accession, so the
 * 1.2 GB source is read but never stored.
 */
process READTHROUGH_ANALYSIS {
    tag "$meta.id"
    label 'process_single'
    maxRetries 2
    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }

    // samtools (reference +4 base) + python, as used by the panel CN check
    conda "${projectDir}/modules/local/cn_candidate_check/environment.yml"

    input:
    tuple val(meta), path(vep_vcf), path(nmd_json), path(domain_json)
    tuple val(ref_meta), path(fasta), path(fai)
    tuple path(clinvar_vcf), path(clinvar_tbi)

    output:
    path "${meta.id}.readthrough.json", emit: json
    path "${meta.id}.readthrough.tsv" , emit: tsv
    path "versions.yml"               , emit: versions

    script:
    // kept out of nextflow.config on purpose: adding params there changed the hash of unrelated
    // reference-prep tasks and forced an expensive rebuild on -resume
    def uniprot = params.containsKey('mechanism_uniprot') ? params.mechanism_uniprot : 'O60566'   // BubR1
    def am_url  = params.containsKey('alphamissense_aa_url') ? params.alphamissense_aa_url :
        'https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz'
    def script_md5 = file("${projectDir}/bin/readthrough_analysis.py").text.md5()
    """
    readthrough_analysis.py --vcf ${vep_vcf} --fasta ${fasta} --nmd ${nmd_json} --domains ${domain_json} \\
        --gene ${params.mechanism_gene} --uniprot ${uniprot} \\
        --am-url '${am_url}' --clinvar ${clinvar_vcf} \\
        --prefix ${meta.id}.readthrough

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
        alphamissense_substitutions: "${am_url}"
        uniprot: "${uniprot}"
        clinvar: "${params.clinvar_release}"
        readthrough_analysis_py_md5: "${script_md5}"
    END_VERSIONS
    """

    stub:
    """
    echo '{"gene": "${params.mechanism_gene}", "ptcs": []}' > ${meta.id}.readthrough.json
    printf 'gene\\thgvsp\\tstop_codon\\n' > ${meta.id}.readthrough.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.11
    END_VERSIONS
    """
}
