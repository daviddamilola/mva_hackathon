/*
 * OptiType 1.5.0 class I typing (HLA-A/B/C), DNA mode, CBC ILP solver, bundled reference.
 * Independent of T1K (different reference and algorithm); used to cross-check class I.
 */
process OPTITYPE_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // pinned env incl. solver; docker/singularity build from environment.yml via Wave

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.optitype_result.tsv"), emit: result
    path "${meta.id}.optitype_coverage.pdf", optional: true  , emit: plot
    path "versions.yml"                                      , emit: versions

    script:
    """
    optitype run -i ${r1} -i ${r2} --dna -o opti_out -p ${meta.id} --solver cbc --threads ${task.cpus}
    cp opti_out/${meta.id}_result.tsv ${meta.id}.optitype_result.tsv
    cp opti_out/${meta.id}_coverage_plot.pdf ${meta.id}.optitype_coverage.pdf 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        optitype: \$(optitype --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
        razers3: "3.5.12"
        solver: "cbc 2.10.13"
    END_VERSIONS
    """

    stub:
    """
    printf '\\tA1\\tA2\\tB1\\tB2\\tC1\\tC2\\tReads\\tObjective\\n0\\tA*01:01\\tA*02:01\\tB*07:02\\tB*08:01\\tC*07:01\\tC*07:02\\t1\\t1\\n' > ${meta.id}.optitype_result.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        optitype: 1.5.0
    END_VERSIONS
    """
}
