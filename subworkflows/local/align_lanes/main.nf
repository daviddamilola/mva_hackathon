/*
 * FASTQ lanes on Hugging Face -> per-lane CRAM (FASTQ deleted after each lane) -> merged, duplicate-marked CRAM.
 */
include { FETCH_ALIGN_LANE } from '../../../modules/local/fetch_align_lane'
include { MERGE_MARKDUP    } from '../../../modules/local/merge_markdup'

workflow ALIGN_LANES {
    take:
    ch_lanes    // [ meta, r1, r1_sha256, r2, r2_sha256 ]
    ch_fasta    // [ meta, fa ]
    ch_fai      // [ meta, fai ]
    ch_bwa      // [ meta, bwa_dir ]

    main:
    ch_ref = ch_fasta
        .join(ch_fai)
        .join(ch_bwa)
        .first()

    FETCH_ALIGN_LANE(ch_lanes, ch_ref)

    ch_lane_crams = FETCH_ALIGN_LANE.out.cram
        .map { meta, cram -> [ [ id: meta.id ], cram ] }
        .groupTuple()

    MERGE_MARKDUP(ch_lane_crams, ch_fasta.join(ch_fai).first())

    emit:
    cram      = MERGE_MARKDUP.out.cram           // [ meta, cram, crai ]
    markdup   = MERGE_MARKDUP.out.stats
    checksums = FETCH_ALIGN_LANE.out.checksums
    versions  = FETCH_ALIGN_LANE.out.versions.first().mix(MERGE_MARKDUP.out.versions)
}
