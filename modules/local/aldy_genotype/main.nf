/*
 * Aldy structural-aware star-allele calling, one gene per task, directly on the CRAM with the
 * explicit reference. Aldy is licensed free for academic/non-commercial use only.
 * A gene Aldy cannot call yields an empty .aldy (reported as no_call), not a pipeline failure.
 */
process ALDY_GENOTYPE {
    tag "${meta.id}:${gene}"
    label 'process_single'
    memory { 4.GB * task.attempt }

    conda "${moduleDir}/environment.yml"
    // no biocontainer for Aldy; docker/singularity profiles build from environment.yml via Wave

    input:
    tuple val(meta), path(cram), path(crai), path(ref_check), val(gene)
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.${gene}.aldy"), emit: aldy
    path "${meta.id}.${gene}.aldy.log"               , emit: log
    path "versions.yml"                              , emit: versions

    script:
    def profile = params.aldy_profile
    """
    set +e
    aldy genotype -p ${profile} -g ${gene} --genome hg38 -r ${fasta} \\
        -o ${meta.id}.${gene}.aldy -l ${meta.id}.${gene}.aldy.log ${cram} > aldy.stdout 2>&1
    rc=\$?
    set -e
    cat aldy.stdout >> ${meta.id}.${gene}.aldy.log
    if [ \$rc -ne 0 ]; then
        echo "ALDY_FAILED exit=\$rc" >> ${meta.id}.${gene}.aldy.log
        touch ${meta.id}.${gene}.aldy
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aldy: \$(python3 -c 'import importlib.metadata as m; print(m.version("aldy"))')
        ortools: \$(python3 -c 'import importlib.metadata as m; print(m.version("ortools"))')
        aldy_args: "genotype -p ${profile} --genome hg38 -r <reference>"
        aldy_license: "Aldy/IURTC: free for academic/non-commercial use"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.${gene}.aldy ${meta.id}.${gene}.aldy.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aldy: 4.8.3
        ortools: 9.15.6755
        aldy_args: "genotype -p ${params.aldy_profile} --genome hg38 -r <reference>"
        aldy_license: "Aldy/IURTC: free for academic/non-commercial use"
    END_VERSIONS
    """
}
