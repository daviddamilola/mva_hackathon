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
include { GENOME_QC         } from '../subworkflows/local/genome_qc'
include { PROVENANCE        } from '../modules/local/provenance'
include { REPORT            } from '../modules/local/report'
include { REPORT_FIGURES    } from '../modules/local/report_figures'
include { OT_QUERY as OT_LITERATURE   } from '../modules/local/ot_query'
include { OT_QUERY as OT_TRACTABILITY } from '../modules/local/ot_query'
include { OT_NEIGHBOURHOOD    } from '../modules/local/ot_neighbourhood'

workflow MVA {
    main:
    ch_versions = Channel.empty()

    ref_meta = [ id: params.ref_name, url: params.ref_fasta_url, md5: params.ref_fasta_md5 ]
    vcf      = file(params.vcf, checkIfExists: true)
    mask_bed = file(params.ref_mask_bed, checkIfExists: true)

    PREPARE_REFERENCE(ref_meta, mask_bed, vcf, !params.skip_alignment)
    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    def sample_meta = [ id: params.sample_id, sex: params.sample_sex ]

    // duplicate-metrics file for QC. When the CRAM is supplied rather than produced here, the file the
    // alignment step published is used, so skipping alignment does not silently drop the duplicate rate
    ch_markdup = Channel.empty()
    if (params.cram) {
        // pre-existing CRAM (e.g. from an earlier run); VALIDATE_CRAM_REF still checks it
        ch_cram = Channel.of([ sample_meta,
                               file(params.cram, checkIfExists: true),
                               file(params.crai ?: "${params.cram}.crai", checkIfExists: true) ])
        def md = file(params.markdup_stats ?: "${params.outdir}/alignment/${params.sample_id}.markdup.stats.txt")
        if (md.exists()) {
            ch_markdup = Channel.value(md)
        } else if (!params.skip_qc) {
            error "Supplying --cram needs the duplicate-metrics file for QC: expected ${md}, or pass --markdup_stats"
        }
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
        ch_markdup = ALIGN_LANES.out.markdup.first()
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

    if (!params.skip_qc && !params.skip_aneuploidy) {
        GENOME_QC(
            PGX_CALL.out.valid_cram, ch_vcf_a, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai,
            ch_markdup,
            ANEUPLOIDY.out.mosdepth_bins_summary.first(),
            ANEUPLOIDY.out.summary,
            Channel.empty()
        )
        ch_versions = ch_versions.mix(GENOME_QC.out.versions)
    }

    if (!params.skip_report && !params.skip_pgx && !params.skip_modifier_panel && !params.skip_aneuploidy && !params.skip_qc) {
        // report inputs and figures first: their version records must be in software_versions.yml before PROVENANCE reads it
        def stage1b = [ "pgx_ot.summary.json", "pgx_ot.verdicts.tsv", "pgx_ot.evidence_levels.tsv" ]
            .collect { file("${params.outdir}/pgx/open_targets/${params.sample_id}.${it}") }.findAll { it.exists() }
        // publication counts for the disease terms: a quantitative measure of unmet need
        OT_LITERATURE(Channel.of([ 'literature', params.unmet_need_disease_ids ]))
        // can the causal gene itself be drugged? the answer is the premise of the candidate list
        OT_TRACTABILITY(Channel.of([ 'tractability', params.mechanism_gene ]))
        // what the drugs that exist near this pathway actually do: computed, not asserted
        OT_NEIGHBOURHOOD(file("${projectDir}/assets/mitotic_drug_targets.tsv", checkIfExists: true))
        ch_versions = ch_versions.mix(OT_LITERATURE.out.versions, OT_TRACTABILITY.out.versions, OT_NEIGHBOURHOOD.out.versions)

        ch_report_in = Channel.empty().mix(
            PGX_CALL.out.calls.map { m, t, j -> [ t, j ] }.flatten(),
            PGX_CALL.out.report.map { it[1] },
            PGX_CALL.out.star_compare.map { it[1] },
            PGX_CALL.out.hla_files,
            PGX_CALL.out.forcecall_summary,
            PGX_CALL.out.mt_rnr1,
            PGX_CALL.out.counterfactual,
            PGX_CALL.out.forcecall.map { m, vcf, csi -> vcf },
            PGX_CALL.out.locus_depth.map { it[1] },
            MODIFIER_PANEL.out.summary,
            MODIFIER_PANEL.out.rare_functional.map { it[1] },
            MODIFIER_PANEL.out.nmd,
            MODIFIER_PANEL.out.domain_map,
            MODIFIER_PANEL.out.domain_figure,
            MODIFIER_PANEL.out.readthrough,
            MODIFIER_PANEL.out.readthrough_tsv,
            ANEUPLOIDY.out.summary.map { it[1] },
            ANEUPLOIDY.out.arms_json,
            ANEUPLOIDY.out.figures.flatten(),
            ANEUPLOIDY.out.cnv_summary,
            ANEUPLOIDY.out.cn_candidates,
            ANEUPLOIDY.out.gene_cn.map { it[1] },
            ANEUPLOIDY.out.bins,
            ANEUPLOIDY.out.exon_cn,
            GENOME_QC.out.qc_json.map { it[1] },
            GENOME_QC.out.samtools_stats,
            GENOME_QC.out.somalier_html,
            GENOME_QC.out.roh_txt,
            OT_LITERATURE.out.json.map { it[2] },
            OT_TRACTABILITY.out.json.map { it[2] },
            OT_NEIGHBOURHOOD.out.json,
            Channel.fromList(stage1b)
        ).collect().map { files -> [ sample_meta, files ] }
        // figures are drawn from the same staged inputs, then handed to the renderer alongside them
        REPORT_FIGURES(ch_report_in, file("${projectDir}/assets/pgx_locus_genes.tsv", checkIfExists: true))
        ch_report_all = ch_report_in.join(REPORT_FIGURES.out.figures.mix(REPORT_FIGURES.out.manifest.map { [ sample_meta, it ] }).groupTuple())
            .map { meta, files, figs -> [ meta, files + figs.flatten() ] }
        ch_versions = ch_versions.mix(REPORT_FIGURES.out.versions)
    }

    ch_sw_versions = ch_versions
        .map { f -> f.text.stripIndent().trim() + '\n' }
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)

    if (!params.skip_report && !params.skip_pgx && !params.skip_modifier_panel && !params.skip_aneuploidy && !params.skip_qc) {
        def git = { List args -> try { return ((['git', '-C', projectDir.toString()] + args) as List<String>).execute().text.trim() } catch (e) { return 'unavailable' } }
        def run_meta = [
            pipeline : [ name: workflow.manifest.name, version: workflow.manifest.version,
                         git_commit: git(['rev-parse', 'HEAD']), git_dirty: git(['status', '--porcelain']) != '',
                         nextflow_version: nextflow.version.toString(), session_id: workflow.sessionId.toString(),
                         profile: workflow.profile, run_start: workflow.start.toString(), command: workflow.commandLine ],
            reference: [ build: params.genome, name: params.ref_name, source: params.ref_fasta_url, source_md5: params.ref_fasta_md5,
                         mask: 'GRC exclusions v2 (inferred), assets/grch38_grc_exclusions_v2_inferred.bed' ],
            databases: [ vep_cache: "Ensembl ${params.vep_cache_version} indexed GRCh38 (offline; checksum-verified)",
                         gnomad: "v${params.gnomad_version} joint (panel regions via ${params.gnomad_base_url})",
                         clinvar: "GRCh38 ${params.clinvar_release} (md5-verified)",
                         alphamissense: "hg38 (${params.alphamissense_url}; md5 ${params.alphamissense_md5_b64})",
                         ipd_imgt_hla: "${params.imgt_hla_version} (${params.imgt_hla_dat_url})",
                         spliceai_annotation: 'Ensembl 116 MANE Select (assets/spliceai_panel_mane_ensembl116.txt)',
                         panel: 'Ensembl 116, +/-10 kb (assets/mechanistic_panel_grch38.bed)',
                         mappability: "Umap k100 (${params.umap_url})", blacklist: "ENCODE v2 (${params.blacklist_url})",
                         segdups: "UCSC genomicSuperDups (${params.segdup_url})", cytoband: "UCSC (${params.cytoband_url})",
                         ancestry_reference: '1000 Genomes phase 3 (somalier profiles)',
                         contamination_resource: params.verifybamid2_resource ],
            params   : params.findAll { k, v -> !(v instanceof Map) }
        ]
        ch_run_meta = Channel.of(groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(run_meta))).collectFile(name: 'run_meta.json')
        def ot_sum = file("${params.outdir}/pgx/open_targets/${params.sample_id}.pgx_ot.summary.json")
        PROVENANCE(ch_sw_versions, ch_run_meta, PGX_CALL.out.report.map { it[1] },
                   ot_sum.exists() ? ot_sum : file("${projectDir}/assets/NO_OT_SUMMARY"))

        REPORT(ch_report_all, PROVENANCE.out.json, file("${projectDir}/assets/report", checkIfExists: true))
    }

    emit:
    versions = ch_versions
}
