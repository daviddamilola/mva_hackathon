/*
 * PyPGx NGS genotyping for one gene (runs in parallel with Aldy). Depth/control inputs are passed
 * for every gene; PyPGx ignores them for genes without SV-defined star alleles.
 * A gene PyPGx cannot call yields a row with an empty genotype (reported), not a pipeline failure.
 */
process PYPGX_GENOTYPE {
    tag "${meta.id}:${gene}"
    label 'process_single'
    memory { 4.GB * task.attempt }

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(vcf), path(tbi), path(depth), path(control), val(gene)
    path bundle

    output:
    tuple val(meta), path("${meta.id}.${gene}.pypgx.tsv"), emit: calls
    path "${meta.id}.${gene}.pypgx.log"                  , emit: log
    path "${meta.id}.${gene}.pypgx_results.zip", optional: true, emit: results
    path "versions.yml"                                  , emit: versions

    script:
    """
    export PYPGX_BUNDLE=\$PWD/${bundle}
    set +e
    pypgx run-ngs-pipeline ${gene} out --variants ${vcf} --depth-of-coverage ${depth} \\
        --control-statistics ${control} --assembly GRCh38 --platform WGS \\
        --do-not-plot-copy-number --do-not-plot-allele-fraction > ${meta.id}.${gene}.pypgx.log 2>&1
    rc=\$?
    set -e
    printf 'Gene\\tGenotype\\tPhenotype\\tHaplotype1\\tHaplotype2\\tCNV\\tStatus\\n' > ${meta.id}.${gene}.pypgx.tsv
    if [ \$rc -eq 0 ] && [ -f out/results.zip ]; then
        cp out/results.zip ${meta.id}.${gene}.pypgx_results.zip
        cnv=\$( [ -f out/cnv-calls.zip ] && pypgx print-data out/cnv-calls.zip | awk -F'\\t' 'NR>1{print \$2}' || echo "" )
        pypgx print-data out/results.zip | awk -F'\\t' -v g=${gene} -v c="\$cnv" 'NR>1{print g"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"c"\\tcalled"}' >> ${meta.id}.${gene}.pypgx.tsv
    else
        printf '%s\\t\\t\\t\\t\\t\\tfailed (exit %s; see log)\\n' ${gene} \$rc >> ${meta.id}.${gene}.pypgx.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: \$(pypgx --version 2>&1 | tail -1 | sed 's/pypgx //')
        pypgx_args: "run-ngs-pipeline --assembly GRCh38 --platform WGS"
    END_VERSIONS
    """

    stub:
    """
    printf 'Gene\\tGenotype\\tPhenotype\\tHaplotype1\\tHaplotype2\\tCNV\\tStatus\\n${gene}\\t*1/*1\\tstub\\t\\t\\t\\tcalled\\n' > ${meta.id}.${gene}.pypgx.tsv
    touch ${meta.id}.${gene}.pypgx.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: 0.27.0
        pypgx_args: "run-ngs-pipeline --assembly GRCh38 --platform WGS"
    END_VERSIONS
    """
}
