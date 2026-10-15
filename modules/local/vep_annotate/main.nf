/*
 * Ensembl VEP, offline with the indexed cache (no network lookups), explicit reference, MANE Select
 * flagged. External annotations are joined as exact-allele custom VCFs: gnomAD v4.1 joint (panel
 * slice), ClinVar (dated release), AlphaMissense (panel slice), SpliceAI (scored locally).
 */
process VEP_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:116.2--pl5321h2a3209d_0' :
        'quay.io/biocontainers/ensembl-vep:116.2--pl5321h2a3209d_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(spliceai), path(spliceai_tbi)
    tuple val(ref_meta), path(fasta), path(fai)
    path vep_cache
    tuple path(gnomad), path(gnomad_tbi)
    tuple path(clinvar), path(clinvar_tbi)
    tuple path(am), path(am_tbi)

    output:
    tuple val(meta), path("${meta.id}.panel.vep.vcf.gz"), emit: vcf
    path "${meta.id}.panel.vep_summary.html"             , emit: html
    path "versions.yml"                                  , emit: versions

    script:
    def gfields = params.gnomad_fields.tokenize(',').collect { it.trim() }.join('%')
    """
    vep --offline --cache --dir_cache ${vep_cache} --cache_version ${params.vep_cache_version} \\
        --species homo_sapiens --assembly GRCh38 --fasta ${fasta} \\
        --input_file ${vcf} --format vcf --vcf --compress_output bgzip \\
        --output_file ${meta.id}.panel.vep.vcf.gz --stats_file ${meta.id}.panel.vep_summary.html \\
        --fork ${task.cpus} --mane --mane_select --canonical --symbol --hgvs --protein --numbers \\
        --variant_class --biotype --check_existing --no_escape \\
        --custom file=${gnomad},short_name=gnomAD_joint,format=vcf,type=exact,coords=0,fields=${gfields} \\
        --custom file=${clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN \\
        --custom file=${am},short_name=AlphaMissense,format=vcf,type=exact,coords=0,fields=AM_PATHOGENICITY%AM_CLASS \\
        --custom file=${spliceai},short_name=SpliceAI,format=vcf,type=exact,coords=0,fields=SpliceAI

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensembl_vep: \$(vep --help 2>&1 | grep -m1 'ensembl-vep' | sed -E 's/.*: *//')
        vep_cache: "homo_sapiens ${params.vep_cache_version}_GRCh38 (indexed, offline)"
        clinvar_release: "${params.clinvar_release}"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.panel.vep.vcf.gz ${meta.id}.panel.vep_summary.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensembl_vep: 116.2
        vep_cache: "homo_sapiens ${params.vep_cache_version}_GRCh38 (indexed, offline)"
        clinvar_release: "${params.clinvar_release}"
    END_VERSIONS
    """
}
