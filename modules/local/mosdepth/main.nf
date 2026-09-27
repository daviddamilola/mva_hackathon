/*
 * mosdepth mean depth over a BED (genome bins or panel exons), CRAM decoded with the explicit
 * reference, reads with MAPQ < params.mosdepth_mapq ignored. Also produces the coverage summary and
 * cumulative distribution used by genome QC.
 */
process MOSDEPTH {
    tag "${meta.id}:${label}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.14--h87be163_2' :
        'quay.io/biocontainers/mosdepth:0.3.14--h87be163_2' }"

    input:
    tuple val(meta), path(cram), path(crai), path(ref_check), val(label), path(bed)
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    tuple val(meta), val(label), path("${meta.id}.${label}.regions.bed.gz"), emit: regions
    path "${meta.id}.${label}.mosdepth.summary.txt"                      , emit: summary
    path "${meta.id}.${label}.mosdepth.global.dist.txt"                  , emit: dist
    path "versions.yml"                                                  , emit: versions

    script:
    def fast = label == 'bins' ? '--fast-mode' : ''
    """
    mosdepth --threads ${task.cpus} --fasta ${fasta} --by ${bed} --mapq ${params.mosdepth_mapq} \\
        --no-per-base ${fast} ${meta.id}.${label} ${cram}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version | sed 's/mosdepth //')
        mosdepth_args: "--by <${label}> --mapq ${params.mosdepth_mapq} ${fast}"
    END_VERSIONS
    """

    stub:
    """
    printf '1\\t0\\t100000\\t45.0\\n' | gzip -c > ${meta.id}.${label}.regions.bed.gz
    touch ${meta.id}.${label}.mosdepth.summary.txt ${meta.id}.${label}.mosdepth.global.dist.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: 0.3.14
        mosdepth_args: "--by <${label}> --mapq ${params.mosdepth_mapq}"
    END_VERSIONS
    """
}
