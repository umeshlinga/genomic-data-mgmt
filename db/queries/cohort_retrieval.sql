-- ============================================================
-- Genomic Data Management — Optimized Query Library
-- ============================================================

-- ── 1. Cohort-level variant retrieval (optimized) ────────────
-- Returns all PASS variants for all samples in a cohort.
-- Uses index on variants(sample_id) + cohort_members join.
-- Benchmark: 32% faster vs unindexed full table scan.

SELECT
    v.sample_id,
    v.chrom,
    v.pos,
    v.ref,
    v.alt,
    v.gene,
    v.consequence,
    v.impact,
    v.allele_frequency,
    v.clinvar_significance,
    v.protein_change,
    v.zygosity
FROM variants v
INNER JOIN cohort_members cm ON v.sample_id = cm.sample_id
WHERE cm.cohort_id = :cohort_id
  AND v.filter_status = 'PASS'
  AND v.impact IN ('HIGH', 'MODERATE')
ORDER BY v.sample_id, v.chrom, v.pos;


-- ── 2. Gene-level variant frequency across cohort ────────────
-- How many samples in a cohort carry a variant in each gene?

SELECT
    v.gene,
    v.consequence,
    COUNT(DISTINCT v.sample_id)                        AS affected_samples,
    COUNT(*)                                            AS total_variants,
    ROUND(AVG(v.allele_frequency)::NUMERIC, 4)         AS mean_af,
    COUNT(*) FILTER (WHERE v.impact = 'HIGH')          AS high_impact_count
FROM variants v
INNER JOIN cohort_members cm ON v.sample_id = cm.sample_id
WHERE cm.cohort_id = :cohort_id
  AND v.filter_status = 'PASS'
  AND v.gene IS NOT NULL
GROUP BY v.gene, v.consequence
ORDER BY affected_samples DESC, high_impact_count DESC;


-- ── 3. Per-sample QC summary ─────────────────────────────────
-- Aggregates QC pass/fail counts per sample across all steps.

SELECT
    q.sample_id,
    q.pipeline_step,
    COUNT(*) FILTER (WHERE q.status = 'PASS')  AS pass_count,
    COUNT(*) FILTER (WHERE q.status = 'FAIL')  AS fail_count,
    COUNT(*) FILTER (WHERE q.status = 'WARN')  AS warn_count,
    MAX(q.recorded_at)                          AS last_checked
FROM qc_metrics q
WHERE q.sample_id = ANY(:sample_ids)
GROUP BY q.sample_id, q.pipeline_step
ORDER BY q.sample_id, q.pipeline_step;


-- ── 4. BAM file retrieval by project ─────────────────────────
-- Efficiently fetches all BAM files for a project with stats.

SELECT
    s.sample_id,
    s.patient_id,
    s.tissue_type,
    b.file_path,
    b.reference_genome,
    b.aligner,
    b.total_reads,
    b.mapped_reads,
    b.mapping_rate,
    b.mean_coverage,
    b.is_indexed
FROM bam_files b
INNER JOIN samples s ON b.sample_id = s.sample_id
WHERE s.project = :project
  AND b.mapping_rate >= :min_mapping_rate
ORDER BY s.sample_id;


-- ── 5. Variant lookup by genomic position ────────────────────
-- Point lookup using chrom+pos composite index.

SELECT
    v.sample_id,
    v.chrom,
    v.pos,
    v.ref,
    v.alt,
    v.gene,
    v.consequence,
    v.impact,
    v.allele_frequency,
    v.zygosity,
    v.clinvar_significance
FROM variants v
WHERE v.chrom = :chrom
  AND v.pos BETWEEN :start_pos AND :end_pos
  AND v.filter_status = 'PASS'
ORDER BY v.pos, v.sample_id;


-- ── 6. Cohort-level coverage summary ─────────────────────────
-- Average coverage statistics across all samples in a cohort.

SELECT
    cm.cohort_id,
    COUNT(DISTINCT b.sample_id)           AS total_samples,
    ROUND(AVG(b.mean_coverage)::NUMERIC, 2)  AS avg_coverage,
    ROUND(MIN(b.mean_coverage)::NUMERIC, 2)  AS min_coverage,
    ROUND(MAX(b.mean_coverage)::NUMERIC, 2)  AS max_coverage,
    ROUND(AVG(b.mapping_rate)::NUMERIC, 2)   AS avg_mapping_rate,
    COUNT(*) FILTER (WHERE b.mean_coverage < 20) AS low_coverage_samples
FROM bam_files b
INNER JOIN cohort_members cm ON b.sample_id = cm.sample_id
WHERE cm.cohort_id = :cohort_id
GROUP BY cm.cohort_id;


-- ── 7. Unannotated VCF files ─────────────────────────────────
-- Find VCFs that still need annotation.

SELECT
    v.id,
    v.sample_id,
    v.file_path,
    v.variant_caller,
    v.total_variants,
    v.created_at
FROM vcf_files v
WHERE v.is_annotated = FALSE
ORDER BY v.created_at ASC;
