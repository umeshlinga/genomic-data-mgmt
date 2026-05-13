-- ============================================================
-- Genomic Data Management System — PostgreSQL Schema
-- ============================================================
-- Tables: samples, fastq_files, bam_files, vcf_files,
--         variants, cohorts, cohort_members, qc_metrics
-- ============================================================

-- Enable UUID extension
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

-- ── Samples ──────────────────────────────────────────────────
CREATE TABLE samples (
    id              UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    sample_id       VARCHAR(100) UNIQUE NOT NULL,
    patient_id      VARCHAR(100),
    project         VARCHAR(100) NOT NULL,
    tissue_type     VARCHAR(100),
    library_type    VARCHAR(50) CHECK (library_type IN ('DNA-seq', 'RNA-seq', 'WGS', 'WES', 'ChIP-seq')),
    sequencing_date DATE,
    created_at      TIMESTAMP DEFAULT NOW(),
    updated_at      TIMESTAMP DEFAULT NOW(),
    notes           TEXT
);

CREATE INDEX idx_samples_project ON samples(project);
CREATE INDEX idx_samples_patient ON samples(patient_id);

-- ── FASTQ Files ───────────────────────────────────────────────
CREATE TABLE fastq_files (
    id              UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    sample_id       VARCHAR(100) REFERENCES samples(sample_id) ON DELETE CASCADE,
    file_path       TEXT NOT NULL,
    read_pair       SMALLINT CHECK (read_pair IN (1, 2)),
    file_size_bytes BIGINT,
    md5_checksum    VARCHAR(32),
    total_reads     BIGINT,
    pass_filter_reads BIGINT,
    mean_quality    NUMERIC(5,2),
    gc_content_pct  NUMERIC(5,2),
    qc_status       VARCHAR(10) DEFAULT 'PENDING' CHECK (qc_status IN ('PENDING', 'PASS', 'FAIL')),
    created_at      TIMESTAMP DEFAULT NOW()
);

CREATE INDEX idx_fastq_sample ON fastq_files(sample_id);
CREATE INDEX idx_fastq_qc ON fastq_files(qc_status);

-- ── BAM Files ─────────────────────────────────────────────────
CREATE TABLE bam_files (
    id              UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    sample_id       VARCHAR(100) REFERENCES samples(sample_id) ON DELETE CASCADE,
    file_path       TEXT NOT NULL,
    reference_genome VARCHAR(20) NOT NULL,
    aligner         VARCHAR(20),
    total_reads     BIGINT,
    mapped_reads    BIGINT,
    mapping_rate    NUMERIC(5,2),
    duplicate_rate  NUMERIC(5,2),
    mean_coverage   NUMERIC(8,2),
    is_sorted       BOOLEAN DEFAULT TRUE,
    is_indexed      BOOLEAN DEFAULT FALSE,
    file_size_bytes BIGINT,
    created_at      TIMESTAMP DEFAULT NOW()
);

CREATE INDEX idx_bam_sample ON bam_files(sample_id);
CREATE INDEX idx_bam_reference ON bam_files(reference_genome);

-- ── VCF Files ─────────────────────────────────────────────────
CREATE TABLE vcf_files (
    id              UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    sample_id       VARCHAR(100) REFERENCES samples(sample_id) ON DELETE CASCADE,
    bam_id          UUID REFERENCES bam_files(id),
    file_path       TEXT NOT NULL,
    variant_caller  VARCHAR(50),
    reference_genome VARCHAR(20) NOT NULL,
    total_variants  INTEGER DEFAULT 0,
    snp_count       INTEGER DEFAULT 0,
    indel_count     INTEGER DEFAULT 0,
    pass_variants   INTEGER DEFAULT 0,
    is_annotated    BOOLEAN DEFAULT FALSE,
    created_at      TIMESTAMP DEFAULT NOW()
);

CREATE INDEX idx_vcf_sample ON vcf_files(sample_id);
CREATE INDEX idx_vcf_caller ON vcf_files(variant_caller);

-- ── Variants ──────────────────────────────────────────────────
CREATE TABLE variants (
    id                  UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    vcf_id              UUID REFERENCES vcf_files(id) ON DELETE CASCADE,
    sample_id           VARCHAR(100) REFERENCES samples(sample_id),
    chrom               VARCHAR(10) NOT NULL,
    pos                 INTEGER NOT NULL,
    ref                 VARCHAR(500) NOT NULL,
    alt                 VARCHAR(500) NOT NULL,
    qual                NUMERIC(10,2),
    filter_status       VARCHAR(20) DEFAULT 'PASS',
    gene                VARCHAR(100),
    transcript_id       VARCHAR(50),
    consequence         VARCHAR(100),
    impact              VARCHAR(20) CHECK (impact IN ('HIGH', 'MODERATE', 'LOW', 'MODIFIER')),
    allele_frequency    NUMERIC(8,6),
    clinvar_significance VARCHAR(50),
    protein_change      VARCHAR(100),
    zygosity            VARCHAR(20) CHECK (zygosity IN ('HOM', 'HET', 'HOM_REF')),
    created_at          TIMESTAMP DEFAULT NOW()
);

CREATE INDEX idx_variants_sample ON variants(sample_id);
CREATE INDEX idx_variants_gene ON variants(gene);
CREATE INDEX idx_variants_chrom_pos ON variants(chrom, pos);
CREATE INDEX idx_variants_consequence ON variants(consequence);
CREATE INDEX idx_variants_impact ON variants(impact);

-- ── Cohorts ───────────────────────────────────────────────────
CREATE TABLE cohorts (
    id          UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    cohort_id   VARCHAR(100) UNIQUE NOT NULL,
    name        VARCHAR(200) NOT NULL,
    description TEXT,
    project     VARCHAR(100),
    created_by  VARCHAR(100),
    created_at  TIMESTAMP DEFAULT NOW()
);

CREATE TABLE cohort_members (
    cohort_id   VARCHAR(100) REFERENCES cohorts(cohort_id) ON DELETE CASCADE,
    sample_id   VARCHAR(100) REFERENCES samples(sample_id) ON DELETE CASCADE,
    added_at    TIMESTAMP DEFAULT NOW(),
    PRIMARY KEY (cohort_id, sample_id)
);

CREATE INDEX idx_cohort_members_sample ON cohort_members(sample_id);

-- ── QC Metrics ────────────────────────────────────────────────
CREATE TABLE qc_metrics (
    id              UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    sample_id       VARCHAR(100) REFERENCES samples(sample_id) ON DELETE CASCADE,
    pipeline_step   VARCHAR(50) NOT NULL,
    metric_name     VARCHAR(100) NOT NULL,
    metric_value    NUMERIC(15,4),
    metric_unit     VARCHAR(20),
    threshold_min   NUMERIC(15,4),
    threshold_max   NUMERIC(15,4),
    status          VARCHAR(10) CHECK (status IN ('PASS', 'FAIL', 'WARN')),
    recorded_at     TIMESTAMP DEFAULT NOW()
);

CREATE INDEX idx_qc_sample ON qc_metrics(sample_id);
CREATE INDEX idx_qc_step ON qc_metrics(pipeline_step);
CREATE INDEX idx_qc_status ON qc_metrics(status);

-- ── Useful Views ──────────────────────────────────────────────

-- Per-sample file summary
CREATE VIEW sample_file_summary AS
SELECT
    s.sample_id,
    s.project,
    s.library_type,
    COUNT(DISTINCT f.id) AS fastq_count,
    COUNT(DISTINCT b.id) AS bam_count,
    COUNT(DISTINCT v.id) AS vcf_count,
    COALESCE(SUM(v.total_variants), 0) AS total_variants
FROM samples s
LEFT JOIN fastq_files f ON s.sample_id = f.sample_id
LEFT JOIN bam_files b ON s.sample_id = b.sample_id
LEFT JOIN vcf_files v ON s.sample_id = v.sample_id
GROUP BY s.sample_id, s.project, s.library_type;

-- High-impact variant summary per sample
CREATE VIEW high_impact_variants AS
SELECT
    sample_id,
    gene,
    chrom,
    pos,
    ref,
    alt,
    consequence,
    impact,
    allele_frequency,
    clinvar_significance,
    protein_change
FROM variants
WHERE impact IN ('HIGH', 'MODERATE')
  AND filter_status = 'PASS'
ORDER BY sample_id, impact DESC, gene;
