/*
 * ClinVar GRCh38 VCF, pinned to a dated weekly release and verified against NCBI's .md5.
 * ClinVar uses no-chr contig names, matching the pipeline reference. storeDir: fetched once.
 */
process CLINVAR_FETCH {
    tag "clinvar_${release}"
    label 'process_single'
    storeDir "${params.reference_dir}/clinvar"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    val release

    output:
    tuple path("clinvar_${release}.vcf.gz"), path("clinvar_${release}.vcf.gz.tbi"), emit: vcf

    script:
    def base = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38'
    """
    for i in 1 2 3 4 5; do
        curl -fsSL -C - --retry 5 -o clinvar_${release}.vcf.gz ${base}/clinvar_${release}.vcf.gz && break || sleep 30
    done
    expected=\$(curl -fsSL ${base}/clinvar_${release}.vcf.gz.md5 | cut -d' ' -f1)
    verify_checksum.py md5 \$expected clinvar_${release}.vcf.gz
    curl -fsSL -o clinvar_${release}.vcf.gz.tbi ${base}/clinvar_${release}.vcf.gz.tbi || tabix -p vcf clinvar_${release}.vcf.gz
    """

    stub:
    """
    touch clinvar_${release}.vcf.gz clinvar_${release}.vcf.gz.tbi
    """
}
