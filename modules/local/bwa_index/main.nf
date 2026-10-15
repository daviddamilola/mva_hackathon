process BWA_INDEX {
    tag "$meta.id"
    label 'process_single'
    label 'process_long'
    memory { 8.GB * task.attempt }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.19--h577a1d6_1' :
        'quay.io/biocontainers/bwa:0.7.19--h577a1d6_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwa"), emit: index
    path "versions.yml"         , emit: versions

    script:
    """
    mkdir bwa
    bwa index -a bwtsw -p bwa/${fasta.baseName} ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | sed -n 's/^Version: //p')
    END_VERSIONS
    """

    stub:
    """
    mkdir bwa
    for ext in amb ann bwt pac sa; do touch bwa/${fasta.baseName}.\$ext; done
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: 0.7.19
    END_VERSIONS
    """
}
