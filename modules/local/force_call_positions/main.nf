/*
 * Force-call arbitrary positions (Open Targets reverse panel) from the CRAM. Unconstrained
 * 'bcftools call -m' writes every targeted position, hom-ref included, with the observed alleles.
 * Nucleotide genotypes are compared downstream, so no allele constraint is needed.
 * An empty site list produces a header-only VCF.
 */
process FORCE_CALL_POSITIONS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.22--h3a4d415_1' :
        'quay.io/biocontainers/bcftools:1.22--h3a4d415_1' }"

    input:
    tuple val(meta), path(cram), path(crai), path(ref_check), path(sites), path(regions)
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.ot_forcecall.vcf.gz"), path("${meta.id}.ot_forcecall.vcf.gz.csi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def sex = meta.sex ?: 'F'
    """
    set -o pipefail
    printf '%s\\t%s\\n' '${meta.id}' '${sex}' > sample_sex.tsv
    if [ -s ${sites} ]; then
        bgzip -c ${sites} > sites.tsv.gz && tabix -s1 -b2 -e2 sites.tsv.gz
        sort -k1,1V -k2,2n ${regions} > regions.sorted.bed
        bcftools mpileup --threads ${task.cpus} --fasta-ref ${fasta} -R regions.sorted.bed -T sites.tsv.gz \\
            -a FORMAT/AD,FORMAT/DP -Ou ${cram} \\
            | bcftools call -m --ploidy GRCh38 -S sample_sex.tsv -Oz -o ${meta.id}.ot_forcecall.vcf.gz
    else
        printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${meta.id}\\n' | bgzip -c > ${meta.id}.ot_forcecall.vcf.gz
    fi
    bcftools index ${meta.id}.ot_forcecall.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | cut -d' ' -f2)
        call_args: "-m --ploidy GRCh38 -S sex (unconstrained)"
    END_VERSIONS
    """

    stub:
    """
    printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${meta.id}\\n' | gzip -c > ${meta.id}.ot_forcecall.vcf.gz
    touch ${meta.id}.ot_forcecall.vcf.gz.csi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.22
        call_args: "-m --ploidy GRCh38 -S sex (unconstrained)"
    END_VERSIONS
    """
}
