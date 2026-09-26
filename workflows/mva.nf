/*
 * Main workflow. Branches implemented so far: reference preparation, alignment, stage 1a PGx calling.
 * Stage 1b (Open Targets annotation) is a separate entry point: pgx_annotate.nf.
 * Planned: modifier panel, aneuploidy burden, genome QC, report (see docs/).
 */
include { PREPARE_REFERENCE } from '../subworkflows/local/prepare_reference'
include { ALIGN_LANES       } from '../subworkflows/local/align_lanes'
include { PGX_CALL          } from '../subworkflows/local/pgx_call'

workflow MVA {
    main:
    ch_versions = Channel.empty()

    ref_meta = [ id: params.ref_name, url: params.ref_fasta_url, md5: params.ref_fasta_md5 ]
    vcf      = file(params.vcf, checkIfExists: true)
    mask_bed = file(params.ref_mask_bed, checkIfExists: true)

    PREPARE_REFERENCE(ref_meta, mask_bed, vcf, !params.skip_alignment)
    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    def sample_meta = [ id: params.sample_id, sex: params.sample_sex ]

    if (params.cram) {
        // pre-existing CRAM (e.g. from an earlier run); VALIDATE_CRAM_REF still checks it
        ch_cram = Channel.of([ sample_meta,
                               file(params.cram, checkIfExists: true),
                               file(params.crai ?: "${params.cram}.crai", checkIfExists: true) ])
    } else if (!params.skip_alignment) {
        ch_lanes = Channel
            .fromPath(params.lanes, checkIfExists: true)
            .splitCsv(header: true, strip: true)
            .map { row ->
                def meta = [ id: params.sample_id, lane: row.lane, flowcell: row.flowcell, sample_index: row.sample_index ]
                [ meta, row.r1, row.r1_sha256, row.r2, row.r2_sha256 ]
            }

        ALIGN_LANES(
            ch_lanes,
            PREPARE_REFERENCE.out.fasta,
            PREPARE_REFERENCE.out.fai,
            PREPARE_REFERENCE.out.bwa
        )
        ch_versions = ch_versions.mix(ALIGN_LANES.out.versions)
        ch_cram = ALIGN_LANES.out.cram.map { meta, cram, crai -> [ sample_meta, cram, crai ] }
    } else {
        ch_cram = Channel.empty()
    }

    if (!params.skip_pgx) {
        PGX_CALL(ch_cram, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai, PREPARE_REFERENCE.out.dict)
        ch_versions = ch_versions.mix(PGX_CALL.out.versions)
    }

    ch_versions
        .map { f -> f.text.stripIndent().trim() + '\n' }
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)

    emit:
    versions = ch_versions
}
