/*
 * Report: Markdown + LaTeX rendered from pipeline outputs with Jinja2 (bin/render_report.py,
 * assets/report). All inputs are staged; the renderer resolves them by file name, and anything missing
 * is printed as missing. The PDF is compiled only if a LaTeX engine is on PATH (report.tex is the
 * deliverable either way).
 */
process REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(inputs, stageAs: 'res/*')
    path provenance
    path templates

    output:
    tuple val(meta), path("${meta.id}.track2_report.md"), path("${meta.id}.track2_report.tex"), emit: report
    path "figures"                                                                          , emit: figures
    path "${meta.id}.track2_report.pdf", optional: true                                      , emit: pdf
    path "versions.yml"                                                                     , emit: versions

    script:
    """
    render_report.py --results res --sample ${meta.id} --templates ${templates} --provenance ${provenance} \\
        --out-md ${meta.id}.track2_report.md --out-tex ${meta.id}.track2_report.tex --figdir figures
    mkdir -p figures
    if command -v latexmk >/dev/null 2>&1; then
        latexmk -pdf -interaction=nonstopmode -halt-on-error ${meta.id}.track2_report.tex > latex.log 2>&1 || echo "LaTeX compile failed (see latex.log); .tex remains the deliverable"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jinja2: \$(python3 -c 'import jinja2; print(jinja2.__version__)')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p figures; touch ${meta.id}.track2_report.md ${meta.id}.track2_report.tex
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jinja2: 3.1.6
    END_VERSIONS
    """
}
