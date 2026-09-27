#!/usr/bin/env nextflow
/*
 * Stage 1b: Open Targets pharmacogenomic annotation of the stage 1a calls, driven by a candidate
 * drug list. Separate entry point so it can be re-run as the Track 2 candidate set changes,
 * without re-running stage 1a:
 *
 *   CONDA_SUBDIR=osx-arm64 nextflow run pgx_annotate.nf -profile conda --drug_list_file my_drugs.yml
 *
 * Drug list: a YAML file naming the candidates, their axis (rationale group) and tier. Names are
 * resolved to ChEMBL IDs through the Open Targets search API at run time (RESOLVE_CANDIDATES), which
 * also applies the anti-target rules; an unresolved name is reported, not dropped. See
 * assets/drugs_track2.template.yml for the schema. The legacy form (a params file carrying
 * 'drugs: [CHEMBL185, ...]') is still accepted and skips name resolution.
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
include { RESOLVE_CANDIDATES    } from './modules/local/resolve_candidates'
include { OT_SITES             } from './modules/local/ot_sites'
include { FORCE_CALL_POSITIONS } from './modules/local/force_call_positions'
include { PGX_JOIN             } from './modules/local/pgx_join'

workflow {
    if (!params.drug_list_file && !params.drugs) {
        error "Provide a candidate drug list: --drug_list_file <yaml> (see assets/drugs_track2.template.yml), " +
              "or the legacy form -params-file with 'drugs: [CHEMBL185, ...]'"
    }

    // the drug list names the candidates; label and example flag are read from the same file, so the
    // report can never show one list's verdicts under another list's name
    def drug_list = params.drug_list_file ? file(params.drug_list_file, checkIfExists: true) : null
    def spec      = drug_list ? new org.yaml.snakeyaml.Yaml().load(drug_list.text) : [:]
    def list_label   = spec.drug_list_label   ?: params.drug_list_label
    def list_example = spec.drug_list_example != null ? spec.drug_list_example : params.drug_list_example
    if (!list_label) {
        error "Name the drug list with drug_list_label (e.g. track2_candidates_v1); it is recorded in the summary and shown in the report"
    }
    if (drug_list && !(spec.candidates)) {
        error "${drug_list}: 'candidates' is empty. Fill it in (see assets/drugs_track2.template.yml); nothing is assumed."
    }
    if (drug_list) {
        def missing = (spec.candidates ?: []).findAll { !it.name || it.axis == null || it.tier == null }
        if (missing) error "${drug_list}: every candidate needs name, axis and tier; incomplete: ${missing}"
        def axes = (spec.axes ?: [:]).keySet().collect { it.toString() }
        def orphan = (spec.candidates ?: []).collect { it.axis.toString() }.unique().findAll { !(it in axes) }
        if (orphan) error "${drug_list}: candidates reference axes that are not defined: ${orphan.join(', ')}"
        log.info "Stage 1b: ${spec.candidates.size()} candidate(s) by name across ${axes.size()} axes, list '${list_label}'"
    }

    def drugs = params.drugs ? (params.drugs instanceof List ? params.drugs : params.drugs.toString().tokenize(','))
        .collect { it.toString().trim().toUpperCase() }.unique() : []
    drugs.each { if (!(it ==~ /CHEMBL\d+/)) error "Not a ChEMBL ID: ${it}" }
    if (drugs) log.info "Stage 1b: ${drugs.size()} drug(s) given as ChEMBL IDs: ${drugs.join(', ')}"

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

    // names -> ChEMBL IDs (Open Targets search) + anti-target filter; only candidates that resolve
    // and survive the filter are queried, and the rest stay in the tables with their status
    if (drug_list) {
        RESOLVE_CANDIDATES(Channel.value([ sample_meta, drug_list ]))
        ch_candidates = RESOLVE_CANDIDATES.out.json.map { m, j -> j }
        ch_query_ids  = RESOLVE_CANDIDATES.out.query_ids
            .splitText().map { it.trim() }.filter { it }
        ch_resolve_versions = RESOLVE_CANDIDATES.out.versions
    } else {
        ch_candidates = Channel.value(file("${projectDir}/assets/no_candidates.json", checkIfExists: true))
        ch_query_ids  = Channel.fromList(drugs)
        ch_resolve_versions = Channel.empty()
    }

    OT_DRUG(ch_query_ids.map { [ 'drug', it ] })
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
        OT_META.out.json.map { it[2] }.first(),
        ch_candidates.first(),
        list_label,
        list_example as boolean
    )

    Channel.empty()
        .mix(VALIDATE_CRAM_REF.out.versions, OT_META.out.versions.first(), OT_SITES.out.versions,
             FORCE_CALL_POSITIONS.out.versions, PGX_JOIN.out.versions, ch_resolve_versions)
        .map { f -> f.text.stripIndent().trim() + '\n' }
        .unique()
        .collectFile(name: 'software_versions_stage1b.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)
}
