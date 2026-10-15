/*
 * Figures for the report, drawn from the pipeline's published outputs (bin/report_figures.py). Every figure
 * is optional: a missing input is recorded in figures.json and the report prints the gap. Nothing is
 * computed here that is not already in a table; the figures show the evidence, they do not create it.
 */
process REPORT_FIGURES {
    tag "$meta.id"
    label 'process_single'

    // numpy + matplotlib, as used by the aneuploidy model; sharing its environment avoids a second identical build
    conda "${projectDir}/modules/local/aneuploidy_model/environment.yml"

    input:
    tuple val(meta), path(inputs, stageAs: 'res/*')
    path locus_genes

    output:
    tuple val(meta), path("figures/*.png"), emit: figures
    path "figures/${meta.id}.figures.json", emit: manifest
    path "versions.yml"                   , emit: versions

    script:
    def script_md5 = file("${projectDir}/bin/report_figures.py").text.md5()
    """
    report_figures.py --results res --sample ${meta.id} --outdir figures --locus-genes ${locus_genes}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)')
        report_figures_py_md5: "${script_md5}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p figures; echo '{}' > figures/${meta.id}.figures.json
    python3 -c 'open("figures/${meta.id}.stub.png","wb").write(b"")'
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.11
    END_VERSIONS
    """
}
