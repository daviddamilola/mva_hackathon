/*
 * Analysis 4: genome context and QC.
 *   SOMALIER (ancestry vs 1000 Genomes, sex stats) | VERIFYBAMID2 (contamination) | BCFTOOLS_ROH (+ het density)
 *   SAMTOOLS_STATS (+ flagstat) | BCFTOOLS_STATS -> QC_SUMMARY (checks with stated thresholds) -> MULTIQC
 */
include { SOMALIER       } from '../../../modules/local/somalier'
include { VERIFYBAMID2   } from '../../../modules/local/verifybamid2'
include { BCFTOOLS_ROH   } from '../../../modules/local/bcftools_roh'
include { SAMTOOLS_STATS } from '../../../modules/local/samtools_stats'
include { BCFTOOLS_STATS } from '../../../modules/local/bcftools_stats'
include { QC_SUMMARY     } from '../../../modules/local/qc_summary'
include { MULTIQC        } from '../../../modules/local/multiqc'

workflow GENOME_QC {
    take:
    ch_valid_cram       // [ meta, cram, crai, ref_check ]
    ch_vcf              // [ meta, vcf, tbi ]
    ch_fasta
    ch_fai
    markdup_stats       // path
    mosdepth_summary    // path (bins)
    aneuploidy_summary  // [ meta, json ]
    extra_multiqc       // channel of paths

    main:
    ch_ref = ch_fasta.join(ch_fai).first()
    SOMALIER(ch_valid_cram, ch_ref)
    VERIFYBAMID2(ch_valid_cram, ch_ref)
    BCFTOOLS_ROH(ch_vcf, ch_ref)
    SAMTOOLS_STATS(ch_valid_cram, ch_ref)
    BCFTOOLS_STATS(ch_vcf)

    ch_qc_in = SAMTOOLS_STATS.out.stats
        .join(SOMALIER.out.results)
        .join(VERIFYBAMID2.out.selfsm)
        .join(BCFTOOLS_ROH.out.summary)
    QC_SUMMARY(ch_qc_in, markdup_stats, mosdepth_summary, aneuploidy_summary.map { it[1] })

    ch_mqc = Channel.empty().mix(
        SAMTOOLS_STATS.out.stats.map { m, s, f -> [ s, f ] }.flatten(),
        BCFTOOLS_STATS.out.stats.map { it[1] },
        VERIFYBAMID2.out.selfsm.map { it[1] },
        SOMALIER.out.results.map { m, a, s -> a },
        mosdepth_summary,
        extra_multiqc
    ).collect()
    MULTIQC(ch_mqc)

    ch_versions = Channel.empty().mix(SOMALIER.out.versions, VERIFYBAMID2.out.versions, BCFTOOLS_ROH.out.versions,
                                      SAMTOOLS_STATS.out.versions, BCFTOOLS_STATS.out.versions, QC_SUMMARY.out.versions,
                                      MULTIQC.out.versions)
    emit:
    qc_json  = QC_SUMMARY.out.json
    multiqc  = MULTIQC.out.report
    versions = ch_versions
}
