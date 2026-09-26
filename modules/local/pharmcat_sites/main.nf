/*
 * PharmCAT's allele-definition positions (shipped with the pinned pharmcat3 package), converted to
 * the reference's no-chr naming, plus the helper files the force-caller needs:
 *   positions.vcf.gz   : -T targets for mpileup (exact positions)
 *   alleles.tsv.gz     : -T for 'bcftools call -C alleles' (CHROM POS REF,ALT)
 *   regions.bed        : merged +/-100 bp windows for -R, so mpileup index-jumps
 *                        instead of streaming the whole CRAM
 */
process PHARMCAT_SITES {
    tag "pharmcat"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/pharmcat3:3.4.0--py313h106432d_0' :
        'quay.io/biocontainers/pharmcat3:3.4.0--py313h106432d_0' }"

    output:
    tuple path("pharmcat_positions.nochr.vcf.gz"), path("pharmcat_positions.nochr.vcf.gz.csi"), path("pharmcat_alleles.tsv.gz"), path("pharmcat_alleles.tsv.gz.tbi"), path("pharmcat_regions.bed"), emit: sites
    path "chr_to_nochr.txt", emit: chr_map
    path "versions.yml"    , emit: versions

    script:
    """
    POS=\$(dirname \$(command -v pharmcat))/pharmcat_positions.vcf.bgz
    for c in \$(seq 1 22) X Y M; do echo "chr\$c \$c"; done > chr_to_nochr.txt

    bcftools annotate --rename-chrs chr_to_nochr.txt \$POS -Oz -o pharmcat_positions.nochr.vcf.gz
    bcftools index pharmcat_positions.nochr.vcf.gz
    bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\n' pharmcat_positions.nochr.vcf.gz | bgzip -c > pharmcat_alleles.tsv.gz
    tabix -s1 -b2 -e2 pharmcat_alleles.tsv.gz
    bcftools query -f '%CHROM\\t%POS0\\t%END\\n' pharmcat_positions.nochr.vcf.gz \\
        | awk 'BEGIN{OFS="\\t"}{s=\$2-100; if(s<0)s=0; print \$1,s,\$3+100}' \\
        | sort -k1,1V -k2,2n \\
        | awk 'BEGIN{OFS="\\t"} \$1==c && \$2<=e {if(\$3>e)e=\$3; next} {if(c!="")print c,s,e; c=\$1;s=\$2;e=\$3} END{print c,s,e}' \\
        > pharmcat_regions.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: \$(pharmcat -version 2>&1 | head -1 | sed 's/PharmCAT //')
        pharmcat_positions: "\$(bcftools view -H pharmcat_positions.nochr.vcf.gz | wc -l | tr -d ' ') positions"
        bcftools: \$(bcftools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    touch pharmcat_positions.nochr.vcf.gz pharmcat_positions.nochr.vcf.gz.csi pharmcat_alleles.tsv.gz pharmcat_alleles.tsv.gz.tbi
    printf '1\\t0\\t200\\n' > pharmcat_regions.bed
    for c in \$(seq 1 22) X Y M; do echo "chr\$c \$c"; done > chr_to_nochr.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: 3.4.0
        pharmcat_positions: "stub"
        bcftools: 1.22
    END_VERSIONS
    """
}
