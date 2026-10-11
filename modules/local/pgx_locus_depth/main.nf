/*
 * Windowed depth over the CYP2D6/CYP2D7 locus (all reads and MAPQ >= 20), the read-level evidence behind
 * the copy-number call that the report's counterfactual turns on. Figure only; no call is made here.
 */
process PGX_LOCUS_DEPTH {
    tag "$meta.id"
    label 'process_single'

    // same tools as the panel CN check (samtools + python); sharing its environment avoids a second identical build
    conda "${projectDir}/modules/local/cn_candidate_check/environment.yml"

    input:
    tuple val(meta), path(cram), path(crai), path(ref_check)
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.pgx_locus_depth.tsv"), emit: depth
    path "versions.yml"                                    , emit: versions

    script:
    // defaults live here rather than in nextflow.config: adding params there changes the hash of
    // unrelated reference-prep tasks and forces an expensive rebuild on -resume
    def region = params.containsKey('pgx_locus_depth_region') ? params.pgx_locus_depth_region
        : '22:42116498-42155810'                                  // CYP2D6 + CYP2D7, GRCh38
    def window = params.containsKey('pgx_locus_depth_window') ? params.pgx_locus_depth_window : 250
    def script_md5 = file("${projectDir}/bin/pgx_locus_depth.py").text.md5()
    """
    pgx_locus_depth.py --cram ${cram} --fasta ${fasta} --region ${region} \\
        --window ${window} --out ${meta.id}.pgx_locus_depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
        python: \$(python3 --version | cut -d' ' -f2)
        pgx_locus_depth_py_md5: "${script_md5}"
    END_VERSIONS
    """

    stub:
    """
    printf 'chrom\\tstart\\tend\\tdepth_all\\tdepth_mq20\\n' > ${meta.id}.pgx_locus_depth.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.22.1
    END_VERSIONS
    """
}
