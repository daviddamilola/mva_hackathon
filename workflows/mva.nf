/*
 * Main workflow. Branches implemented so far: reference preparation and alignment.
 * Planned: PGx, modifier panel, aneuploidy burden, genome QC, report (see docs/).
 */
include { PREPARE_REFERENCE } from '../subworkflows/local/prepare_reference'
include { ALIGN_LANES       } from '../subworkflows/local/align_lanes'

workflow MVA {
    main:
    ch_versions = Channel.empty()

    ref_meta = [ id: params.ref_name, url: params.ref_fasta_url, md5: params.ref_fasta_md5 ]
    vcf      = file(params.vcf, checkIfExists: true)
    mask_bed = file(params.ref_mask_bed, checkIfExists: true)

    PREPARE_REFERENCE(ref_meta, mask_bed, vcf, !params.skip_alignment)
    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    if (!params.skip_alignment) {
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
    }

    ch_versions
        .map { f -> f.text.stripIndent().trim() + '\n' }
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)

    emit:
    versions = ch_versions
}
