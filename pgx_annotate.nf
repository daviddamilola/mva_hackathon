#!/usr/bin/env nextflow
/*
 * Stage 1b: Open Targets pharmacogenomic annotation of the stage 1a calls, driven by a candidate
 * drug list. Separate entry point so it can be re-run as the Track 2 candidate set changes,
 * without re-running stage 1a:
 *
 *   CONDA_SUBDIR=osx-arm64 nextflow run pgx_annotate.nf -profile conda -params-file my_drugs.yml
 *
 * params file:  drugs: [ CHEMBL185, CHEMBL1200751, ... ]
 *
 * Forward: every called gene -> target(ensemblId){pharmacogenomics}, joined to its diplotype.
 * Reverse: every drug -> drug(chemblId){pharmacogenomics}. The union of returned positions is
 *          force-called from the CRAM where not already genotyped, then joined. The drug list
 *          defines the panel; no fixed gene list is used.
 * Output : per-candidate-drug verdict table, forward annotations, evidence-level counts,
 *          join-coverage summary, Open Targets release in provenance.
 * Research hypothesis only; not a treatment recommendation.
 */
nextflow.enable.dsl = 2

include { VALIDATE_CRAM_REF    } from './modules/local/validate_cram_ref'
include { OT_QUERY as OT_META  } from './modules/local/ot_query'
include { OT_QUERY as OT_DRUG  } from './modules/local/ot_query'
include { OT_QUERY as OT_GENES } from './modules/local/ot_query'
include { OT_SITES             } from './modules/local/ot_sites'
include { FORCE_CALL_POSITIONS } from './modules/local/force_call_positions'
include { PGX_JOIN             } from './modules/local/pgx_join'

workflow {
    if (!params.drugs) {
        error "Provide candidate drugs as ChEMBL IDs, e.g. -params-file drugs.yml with 'drugs: [CHEMBL185, ...]'"
    }
    if (!params.drug_list_label) {
        error "Name the drug list with drug_list_label (e.g. track2_candidates_v1); it is recorded in the summary and shown in the report"
    }
    def drugs = (params.drugs instanceof List ? params.drugs : params.drugs.toString().tokenize(','))
        .collect { it.toString().trim().toUpperCase() }.unique()
    drugs.each { if (!(it ==~ /CHEMBL\d+/)) error "Not a ChEMBL ID: ${it}" }
    log.info "Stage 1b: ${drugs.size()} candidate drug(s): ${drugs.join(', ')}"

    def sample_meta = [ id: params.sample_id, sex: params.sample_sex ]
    def pgx_dir     = "${params.outdir}/pgx"
    def ref_prefix  = "${params.reference_dir}/${params.ref_name}"

    ch_cram = Channel.of([ sample_meta,
                           file(params.cram ?: "${params.outdir}/alignment/${params.sample_id}.cram", checkIfExists: true),
                           file(params.crai ?: "${params.cram ?: "${params.outdir}/alignment/${params.sample_id}.cram"}.crai", checkIfExists: true) ])
    ch_ref      = Channel.value([ [ id: params.ref_name ], file("${ref_prefix}.fa", checkIfExists: true), file("${ref_prefix}.fa.fai", checkIfExists: true) ])
    ch_ref_dict = Channel.value([ [ id: params.ref_name ], file("${ref_prefix}.fa", checkIfExists: true),
                                  file("${ref_prefix}.fa.fai", checkIfExists: true), file("${ref_prefix}.dict", checkIfExists: true) ])

    calls_json   = file(params.pgx_calls_json   ?: "${pgx_dir}/${params.sample_id}.pgx_calls.json", checkIfExists: true)
    report_json  = file(params.pharmcat_report  ?: "${pgx_dir}/pharmcat/${params.sample_id}.report.json", checkIfExists: true)
    stage1a_vcf  = file(params.pgx_forcecall_vcf ?: "${pgx_dir}/forcecall/${params.sample_id}.pgx_forcecall.vcf.gz", checkIfExists: true)
    stage1a_csi  = file("${stage1a_vcf}.csi", checkIfExists: true)

    // first node: CRAM <-> reference M5 validation
    VALIDATE_CRAM_REF(ch_cram, ch_ref_dict)

    // called genes come from stage 1a, not a hardcoded list
    def called_genes = new groovy.json.JsonSlurper().parse(calls_json).calls
        .findAll { it.pharmcat_diplotype || it.aldy_diplotype }
        .collect { it.gene }
    log.info "Stage 1b: ${called_genes.size()} called gene(s) from stage 1a"

    OT_META(Channel.of([ 'meta', 'release' ]))
    OT_DRUG(Channel.fromList(drugs).map { [ 'drug', it ] })
    OT_GENES(Channel.of([ 'target', called_genes.join(',') ]))

    ch_ot_jsons = OT_DRUG.out.json.map { it[2] }.mix(OT_GENES.out.json.map { it[2] }).collect()
    OT_SITES(Channel.value([ sample_meta, stage1a_vcf, stage1a_csi ]), ch_ot_jsons)

    FORCE_CALL_POSITIONS(
        VALIDATE_CRAM_REF.out.cram.join(OT_SITES.out.sites),
        ch_ref
    )

    PGX_JOIN(
        FORCE_CALL_POSITIONS.out.vcf.map { meta, vcf, csi -> [ meta, calls_json, report_json, stage1a_vcf, vcf ] },
        OT_DRUG.out.json.map { it[2] }.collect(),
        OT_GENES.out.json.map { it[2] }.first(),
        OT_META.out.json.map { it[2] }.first()
    )

    Channel.empty()
        .mix(VALIDATE_CRAM_REF.out.versions, OT_META.out.versions.first(), OT_SITES.out.versions,
             FORCE_CALL_POSITIONS.out.versions, PGX_JOIN.out.versions)
        .map { f -> f.text.stripIndent().trim() + '\n' }
        .unique()
        .collectFile(name: 'software_versions_stage1b.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)
}
