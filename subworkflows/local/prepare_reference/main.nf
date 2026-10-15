/*
 * Rebuild the exact reference the subject VCF was called against, index it for bwa, and check that
 * its contigs match the VCF header.
 */
include { PREPARE_FASTA       } from '../../../modules/local/prepare_fasta'
include { SAMTOOLS_FAIDX_DICT } from '../../../modules/local/samtools_faidx_dict'
include { BWA_INDEX           } from '../../../modules/local/bwa_index'
include { CHECK_CONTIGS       } from '../../../modules/local/check_contigs'

workflow PREPARE_REFERENCE {
    take:
    ref_meta    // val: [ id, url, md5 ]
    mask_bed    // path
    vcf         // path: subject VCF (header only is read)
    build_bwa   // val: boolean

    main:
    ch_versions = Channel.empty()

    PREPARE_FASTA(ref_meta, mask_bed)
    SAMTOOLS_FAIDX_DICT(PREPARE_FASTA.out.fasta)
    CHECK_CONTIGS(SAMTOOLS_FAIDX_DICT.out.fai, vcf)

    ch_bwa = Channel.empty()
    if (build_bwa) {
        BWA_INDEX(PREPARE_FASTA.out.fasta)
        ch_bwa = BWA_INDEX.out.index
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
    }

    ch_versions = ch_versions.mix(
        PREPARE_FASTA.out.versions,
        SAMTOOLS_FAIDX_DICT.out.versions,
        CHECK_CONTIGS.out.versions
    )

    emit:
    fasta        = PREPARE_FASTA.out.fasta        // [ meta, fa ]
    fai          = SAMTOOLS_FAIDX_DICT.out.fai    // [ meta, fai ]
    dict         = SAMTOOLS_FAIDX_DICT.out.dict   // [ meta, dict ]
    bwa          = ch_bwa                         // [ meta, dir ]
    contig_check = CHECK_CONTIGS.out.report
    versions     = ch_versions
}
