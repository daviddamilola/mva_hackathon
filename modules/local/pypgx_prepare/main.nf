/*
 * Inputs for PyPGx, computed once and shared by every gene:
 *  - slice.bam : PyPGx gene regions cut from the CRAM with the EXPLICIT reference. PyPGx's depth
 *                step has no reference option, so feeding it a CRAM would make htslib resolve the
 *                reference implicitly from the header UR tag.
 *  - input.vcf.gz (create-input-vcf), depth.zip (prepare-depth-of-coverage, restricted to our
 *    genes: PyPGx's default set includes GSTT1 on an ALT contig absent from the no-ALT
 *    reference), control.zip (compute-control-statistics on the control gene).
 */
process PYPGX_PREPARE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // no single biocontainer ships pypgx + pinned pandas + samtools; docker/singularity build via Wave

    input:
    tuple val(meta), path(cram), path(crai), path(ref_check)
    tuple val(ref_meta), path(fasta), path(fai)
    path bundle

    output:
    tuple val(meta), path("input.vcf.gz"), path("input.vcf.gz.tbi"), path("depth.zip"), path("control.zip"), emit: inputs
    path "versions.yml", emit: versions

    script:
    def genes = params.pypgx_genes.tokenize(',').collect { it.trim() }.join(' ')
    """
    set -o pipefail
    export PYPGX_BUNDLE=\$PWD/${bundle}
    pypgx create-regions-bed --assembly GRCh38 | sed 's/^chr//' | sort -k1,1V -k2,2n > regions.bed
    samtools view -@ ${task.cpus} --reference ${fasta} -b -L regions.bed -o slice.bam ${cram}
    samtools index slice.bam

    pypgx create-input-vcf input.vcf.gz ${fasta} slice.bam --assembly GRCh38 --genes ${genes}
    pypgx prepare-depth-of-coverage depth.zip slice.bam --assembly GRCh38 --genes ${genes}
    pypgx compute-control-statistics --assembly GRCh38 ${params.pypgx_control_gene} control.zip slice.bam
    rm slice.bam slice.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: \$(pypgx --version 2>&1 | tail -1 | sed 's/pypgx //')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
        pypgx_bundle: "${bundle}"
        control_gene: "${params.pypgx_control_gene}"
    END_VERSIONS
    """

    stub:
    """
    touch input.vcf.gz input.vcf.gz.tbi depth.zip control.zip
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: 0.27.0
        pandas: 2.3.3
        samtools: 1.22.1
        pypgx_bundle: "${bundle}"
        control_gene: "${params.pypgx_control_gene}"
    END_VERSIONS
    """
}
