/*
 * Analysis 3: mosaic aneuploidy burden (the biomarker), plus panel-gene CNV/SV to close the
 * modifier-panel gap.
 *   FETCH_RESOURCE (Umap k100, ENCODE blacklist v2, UCSC cytoBand, UCSC segdups, Delly hg38 exclude; sha256-pinned)
 *   BIN_ANNOTATE (100 kb bins + panel MANE exons: GC, mappability, blacklist)
 *   MOSDEPTH (bins, fast mode | panel exons)            <- validated CRAM, explicit reference
 *   ANEUPLOIDY_MODEL (GC-corrected LRR + beta-binomial BAF -> per-arm cell fraction, CI, detection limit)
 *   DELLY_CALL (genome-wide SVs) -> GENE_CN (exon CN + SVs overlapping panel genes)
 *   CN_CANDIDATE_CHECK (depth step in the reads for each exon-run CN candidate)
 */
include { FETCH_RESOURCE                } from '../../../modules/local/fetch_resource'
include { BIN_ANNOTATE                  } from '../../../modules/local/bin_annotate'
include { MOSDEPTH as MOSDEPTH_BINS     } from '../../../modules/local/mosdepth'
include { MOSDEPTH as MOSDEPTH_EXONS    } from '../../../modules/local/mosdepth'
include { ANEUPLOIDY_MODEL              } from '../../../modules/local/aneuploidy_model'
include { DELLY_CALL                    } from '../../../modules/local/delly_call'
include { GENE_CN                       } from '../../../modules/local/gene_cn'
include { CN_CANDIDATE_CHECK            } from '../../../modules/local/cn_candidate_check'

workflow ANEUPLOIDY {
    take:
    ch_valid_cram   // [ meta, cram, crai, ref_check ]
    ch_vcf          // [ meta, vcf, tbi ]
    ch_fasta        // [ ref_meta, fa ]
    ch_fai          // [ ref_meta, fai ]

    main:
    ch_ref = ch_fasta.join(ch_fai).first()

    FETCH_RESOURCE(Channel.fromList([
        [ 'umap_k100',         params.umap_url,          params.umap_sha256 ],
        [ 'blacklist_v2',      params.blacklist_url,     params.blacklist_sha256 ],
        [ 'cytoband',          params.cytoband_url,      params.cytoband_sha256 ],
        [ 'delly_exclude',     params.delly_exclude_url, params.delly_exclude_sha256 ],
        [ 'segdup',            params.segdup_url,        params.segdup_sha256 ]
    ]))
    res = FETCH_RESOURCE.out.file.branch { name, f ->
        umap: name == 'umap_k100'; bl: name == 'blacklist_v2'; cyto: name == 'cytoband'; excl: name == 'delly_exclude'
        segdup: name == 'segdup'
    }

    BIN_ANNOTATE(ch_ref, res.umap.first(), res.bl.first(), file(params.spliceai_annotation, checkIfExists: true))

    MOSDEPTH_BINS(ch_valid_cram.combine(BIN_ANNOTATE.out.bins_bed).map { m, c, i, r, b -> [ m, c, i, r, 'bins', b ] }, ch_ref)
    MOSDEPTH_EXONS(ch_valid_cram.combine(BIN_ANNOTATE.out.exons_bed).map { m, c, i, r, b -> [ m, c, i, r, 'panel_exons', b ] }, ch_ref)

    ANEUPLOIDY_MODEL(MOSDEPTH_BINS.out.regions, BIN_ANNOTATE.out.bins, ch_vcf.first(), res.cyto.first(), res.segdup.first())

    DELLY_CALL(ch_valid_cram, ch_ref, res.excl.first())
    GENE_CN(MOSDEPTH_EXONS.out.regions, BIN_ANNOTATE.out.exons, ANEUPLOIDY_MODEL.out.bins_depth, ANEUPLOIDY_MODEL.out.bins,
            DELLY_CALL.out.vcf.first(), file(params.gene_panel_bed, checkIfExists: true), ANEUPLOIDY_MODEL.out.arms_json)
    CN_CANDIDATE_CHECK(ch_valid_cram.join(GENE_CN.out.genes), GENE_CN.out.exons, GENE_CN.out.sv, ch_ref)

    ch_versions = Channel.empty().mix(BIN_ANNOTATE.out.versions, MOSDEPTH_BINS.out.versions, ANEUPLOIDY_MODEL.out.versions,
                                      DELLY_CALL.out.versions, GENE_CN.out.versions,
                                      CN_CANDIDATE_CHECK.out.versions)

    emit:
    summary  = ANEUPLOIDY_MODEL.out.summary
    arms     = ANEUPLOIDY_MODEL.out.arms
    figures  = ANEUPLOIDY_MODEL.out.figures
    gene_cn  = GENE_CN.out.genes
    arms_json = ANEUPLOIDY_MODEL.out.arms_json
    cnv_summary = GENE_CN.out.summary
    cn_candidates = CN_CANDIDATE_CHECK.out.json
    mosdepth_bins_summary = MOSDEPTH_BINS.out.summary
    coverage = MOSDEPTH_BINS.out.summary
    versions = ch_versions
}
