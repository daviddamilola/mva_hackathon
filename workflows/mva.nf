/*
 * Main workflow. Branches implemented so far: reference preparation, alignment, stage 1a PGx calling,
 * modifier panel.
 * Stage 1b (Open Targets annotation) is a separate entry point: pgx_annotate.nf.
 * aneuploidy burden + panel CNV/SV. Planned: genome QC, report (see docs/).
 */
include { PREPARE_REFERENCE } from '../subworkflows/local/prepare_reference'
include { ALIGN_LANES       } from '../subworkflows/local/align_lanes'
include { PGX_CALL          } from '../subworkflows/local/pgx_call'
include { MODIFIER_PANEL    } from '../subworkflows/local/modifier_panel'
include { ANEUPLOIDY        } from '../subworkflows/local/aneuploidy'

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

    if (!params.skip_modifier_panel) {
        ch_vcf = Channel.of([ sample_meta, file(params.vcf, checkIfExists: true), file(params.vcf_tbi ?: "${params.vcf}.tbi", checkIfExists: true) ])
        MODIFIER_PANEL(ch_vcf, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai)
        ch_versions = ch_versions.mix(MODIFIER_PANEL.out.versions)
    }

    if (!params.skip_aneuploidy) {
        if (params.skip_pgx) error "skip_aneuploidy=false needs the validated CRAM from the PGx branch (skip_pgx=false)"
        ch_vcf_a = Channel.of([ sample_meta, file(params.vcf, checkIfExists: true), file(params.vcf_tbi ?: "${params.vcf}.tbi", checkIfExists: true) ])
        ANEUPLOIDY(PGX_CALL.out.valid_cram, ch_vcf_a, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai)
        ch_versions = ch_versions.mix(ANEUPLOIDY.out.versions)
    }

    ch_versions
        .map { f -> f.text.stripIndent().trim() + '\n' }
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)

    emit:
    versions = ch_versions
}
