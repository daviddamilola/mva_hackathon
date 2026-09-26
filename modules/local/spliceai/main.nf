/*
 * SpliceAI run locally on the panel variants only (no external scoring service: subject variants
 * never leave the machine). Explicit reference, built-in GRCh38 annotation, masked scores (-M 1,
 * recommended for variant interpretation), distance -D 50.
 * Model weights: CC BY-NC 4.0; code: GPLv3.
 */
process SPLICEAI {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // no single-image pin with our TF/numpy pins; docker/singularity build from environment.yml via Wave

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.spliceai.vcf.gz"), path("${meta.id}.spliceai.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    export TF_CPP_MIN_LOG_LEVEL=2
    spliceai -I ${vcf} -O spliceai.vcf -R ${fasta} -A ${params.spliceai_annotation} \\
        -D ${params.spliceai_distance} -M ${params.spliceai_mask}
    bgzip -c spliceai.vcf > ${meta.id}.spliceai.vcf.gz
    tabix -p vcf ${meta.id}.spliceai.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spliceai: \$(python3 -c 'import importlib.metadata as m; print(m.version("spliceai"))')
        tensorflow: \$(python3 -c 'import importlib.metadata as m; print(m.version("tensorflow"))')
        spliceai_args: "-A ${params.spliceai_annotation} -D ${params.spliceai_distance} -M ${params.spliceai_mask}"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.spliceai.vcf.gz ${meta.id}.spliceai.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spliceai: 1.3.1
        tensorflow: 2.21.0
        spliceai_args: "-A ${params.spliceai_annotation} -D ${params.spliceai_distance} -M ${params.spliceai_mask}"
    END_VERSIONS
    """
}
