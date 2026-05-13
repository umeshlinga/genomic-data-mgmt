"""
db_manager.py
-------------
Database Manager Module

Handles PostgreSQL connections, FASTQ/BAM/VCF ingestion,
cohort management, and optimized variant retrieval.
"""

import logging
import json
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any
from datetime import date

try:
    import psycopg2
    import psycopg2.extras
    from psycopg2.pool import ThreadedConnectionPool
except ImportError:
    psycopg2 = None

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


@dataclass
class SampleRecord:
    """Represents a genomic sample record."""
    sample_id: str
    project: str
    library_type: str
    patient_id: Optional[str] = None
    tissue_type: Optional[str] = None
    sequencing_date: Optional[date] = None
    notes: Optional[str] = None


@dataclass
class VariantRecord:
    """Represents a single variant for database insertion."""
    sample_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str = ""
    consequence: str = ""
    impact: str = "MODIFIER"
    allele_frequency: float = 0.0
    filter_status: str = "PASS"
    zygosity: str = "HET"
    clinvar_significance: str = ""
    protein_change: str = ""
    vcf_id: Optional[str] = None


class GenomicDBManager:
    """
    Manages all genomic data storage and retrieval operations.

    Provides methods for ingesting FASTQ/BAM/VCF metadata,
    managing cohorts, and running optimized variant queries.

    Parameters
    ----------
    host : str
        PostgreSQL host.
    port : int
        PostgreSQL port (default: 5432).
    database : str
        Database name.
    user : str
        Database user.
    password : str
        Database password.
    pool_size : int
        Connection pool size (default: 5).
    """

    def __init__(
        self,
        host: str,
        port: int = 5432,
        database: str = "genomics",
        user: str = "postgres",
        password: str = "",
        pool_size: int = 5,
    ):
        self.dsn = {
            "host": host,
            "port": port,
            "dbname": database,
            "user": user,
            "password": password,
        }
        self._pool = None
        self._pool_size = pool_size

    def connect(self):
        """Initialize the connection pool."""
        if psycopg2 is None:
            raise ImportError("psycopg2 is required: pip install psycopg2-binary")
        self._pool = ThreadedConnectionPool(
            minconn=1,
            maxconn=self._pool_size,
            **self.dsn,
        )
        logger.info(f"Connected to database: {self.dsn['dbname']}@{self.dsn['host']}")

    def disconnect(self):
        """Close all connections in the pool."""
        if self._pool:
            self._pool.closeall()
            logger.info("Database connections closed")

    @contextmanager
    def _get_connection(self):
        """Context manager for connection pool checkout."""
        conn = self._pool.getconn()
        try:
            yield conn
            conn.commit()
        except Exception:
            conn.rollback()
            raise
        finally:
            self._pool.putconn(conn)

    # ── Sample Management ─────────────────────────────────────

    def upsert_sample(self, sample: SampleRecord) -> str:
        """
        Insert or update a sample record.

        Parameters
        ----------
        sample : SampleRecord

        Returns
        -------
        str
            The sample_id.
        """
        sql = """
            INSERT INTO samples (
                sample_id, project, library_type, patient_id,
                tissue_type, sequencing_date, notes, updated_at
            ) VALUES (
                %(sample_id)s, %(project)s, %(library_type)s, %(patient_id)s,
                %(tissue_type)s, %(sequencing_date)s, %(notes)s, NOW()
            )
            ON CONFLICT (sample_id) DO UPDATE SET
                project = EXCLUDED.project,
                library_type = EXCLUDED.library_type,
                patient_id = EXCLUDED.patient_id,
                tissue_type = EXCLUDED.tissue_type,
                updated_at = NOW()
            RETURNING sample_id;
        """
        with self._get_connection() as conn:
            with conn.cursor() as cur:
                cur.execute(sql, sample.__dict__)
                result = cur.fetchone()
        logger.info(f"Upserted sample: {sample.sample_id}")
        return result[0]

    def get_sample(self, sample_id: str) -> Optional[Dict]:
        """Retrieve a sample record by ID."""
        sql = "SELECT * FROM samples WHERE sample_id = %s;"
        with self._get_connection() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                cur.execute(sql, (sample_id,))
                row = cur.fetchone()
        return dict(row) if row else None

    # ── File Ingestion ────────────────────────────────────────

    def ingest_fastq(
        self,
        sample_id: str,
        file_path: str,
        read_pair: int = 1,
        total_reads: int = 0,
        pass_filter_reads: int = 0,
        mean_quality: float = 0.0,
        gc_content_pct: float = 0.0,
        qc_status: str = "PENDING",
        file_size_bytes: int = 0,
    ) -> str:
        """
        Register a FASTQ file and its QC metrics.

        Returns
        -------
        str
            UUID of the inserted record.
        """
        sql = """
            INSERT INTO fastq_files (
                sample_id, file_path, read_pair, total_reads,
                pass_filter_reads, mean_quality, gc_content_pct,
                qc_status, file_size_bytes
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
            RETURNING id;
        """
        with self._get_connection() as conn:
            with conn.cursor() as cur:
                cur.execute(sql, (
                    sample_id, file_path, read_pair, total_reads,
                    pass_filter_reads, mean_quality, gc_content_pct,
                    qc_status, file_size_bytes,
                ))
                record_id = str(cur.fetchone()[0])
        logger.info(f"Ingested FASTQ: {file_path} → {record_id}")
        return record_id

    def ingest_bam(
        self,
        sample_id: str,
        file_path: str,
        reference_genome: str,
        aligner: str = "bwa",
        total_reads: int = 0,
        mapped_reads: int = 0,
        mapping_rate: float = 0.0,
        mean_coverage: float = 0.0,
        is_indexed: bool = True,
        file_size_bytes: int = 0,
    ) -> str:
        """Register a BAM file and its alignment metrics."""
        sql = """
            INSERT INTO bam_files (
                sample_id, file_path, reference_genome, aligner,
                total_reads, mapped_reads, mapping_rate,
                mean_coverage, is_indexed, file_size_bytes
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            RETURNING id;
        """
        with self._get_connection() as conn:
            with conn.cursor() as cur:
                cur.execute(sql, (
                    sample_id, file_path, reference_genome, aligner,
                    total_reads, mapped_reads, mapping_rate,
                    mean_coverage, is_indexed, file_size_bytes,
                ))
                record_id = str(cur.fetchone()[0])
        logger.info(f"Ingested BAM: {file_path} → {record_id}")
        return record_id

    def ingest_vcf(
        self,
        sample_id: str,
        file_path: str,
        reference_genome: str,
        variant_caller: str = "gatk",
        total_variants: int = 0,
        snp_count: int = 0,
        indel_count: int = 0,
        pass_variants: int = 0,
        bam_id: Optional[str] = None,
    ) -> str:
        """Register a VCF file and its variant statistics."""
        sql = """
            INSERT INTO vcf_files (
                sample_id, file_path, reference_genome, variant_caller,
                total_variants, snp_count, indel_count, pass_variants, bam_id
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
            RETURNING id;
        """
        with self._get_connection() as conn:
            with conn.cursor() as cur:
                cur.execute(sql, (
                    sample_id, file_path, reference_genome, variant_caller,
                    total_variants, snp_count, indel_count, pass_variants, bam_id,
                ))
                record_id = str(cur.fetchone()[0])
        logger.info(f"Ingested VCF: {file_path} → {record_id}")
        return record_id

    # ── Variant Operations ────────────────────────────────────

    def bulk_insert_variants(
        self, variants: List[VariantRecord], batch_size: int = 1000
    ) -> int:
        """
        Bulk insert variants using psycopg2 execute_values for performance.

        Parameters
        ----------
        variants : list of VariantRecord
        batch_size : int
            Number of rows per INSERT batch (default: 1000).

        Returns
        -------
        int
            Total number of variants inserted.
        """
        if not variants:
            return 0

        sql = """
            INSERT INTO variants (
                vcf_id, sample_id, chrom, pos, ref, alt,
                gene, consequence, impact, allele_frequency,
                filter_status, zygosity, clinvar_significance, protein_change
            ) VALUES %s
            ON CONFLICT DO NOTHING;
        """
        rows = [
            (
                v.vcf_id, v.sample_id, v.chrom, v.pos, v.ref, v.alt,
                v.gene, v.consequence, v.impact, v.allele_frequency,
                v.filter_status, v.zygosity, v.clinvar_significance, v.protein_change,
            )
            for v in variants
        ]

        total = 0
        with self._get_connection() as conn:
            with conn.cursor() as cur:
                for i in range(0, len(rows), batch_size):
                    batch = rows[i:i + batch_size]
                    psycopg2.extras.execute_values(cur, sql, batch)
                    total += len(batch)

        logger.info(f"Inserted {total} variants")
        return total

    def get_cohort_variants(
        self,
        cohort_id: str,
        impact_filter: Optional[List[str]] = None,
        gene_filter: Optional[List[str]] = None,
    ) -> List[Dict]:
        """
        Retrieve all PASS variants for a cohort with optional filters.

        Uses optimized index on variants(sample_id) + cohort join.

        Parameters
        ----------
        cohort_id : str
            Cohort identifier.
        impact_filter : list, optional
            Filter by impact levels (e.g. ['HIGH', 'MODERATE']).
        gene_filter : list, optional
            Filter to specific genes.

        Returns
        -------
        list of dict
        """
        conditions = ["cm.cohort_id = %s", "v.filter_status = 'PASS'"]
        params: List[Any] = [cohort_id]

        if impact_filter:
            conditions.append(f"v.impact = ANY(%s)")
            params.append(impact_filter)
        if gene_filter:
            conditions.append(f"v.gene = ANY(%s)")
            params.append(gene_filter)

        where = " AND ".join(conditions)
        sql = f"""
            SELECT
                v.sample_id, v.chrom, v.pos, v.ref, v.alt,
                v.gene, v.consequence, v.impact,
                v.allele_frequency, v.zygosity,
                v.clinvar_significance, v.protein_change
            FROM variants v
            INNER JOIN cohort_members cm ON v.sample_id = cm.sample_id
            WHERE {where}
            ORDER BY v.sample_id, v.chrom, v.pos;
        """
        with self._get_connection() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                cur.execute(sql, params)
                rows = cur.fetchall()
        return [dict(r) for r in rows]

    # ── Cohort Management ─────────────────────────────────────

    def create_cohort(
        self, cohort_id: str, name: str, sample_ids: List[str],
        description: str = "", project: str = ""
    ) -> str:
        """
        Create a new cohort and add samples to it.

        Parameters
        ----------
        cohort_id : str
        name : str
        sample_ids : list of str
        description : str
        project : str

        Returns
        -------
        str
            The cohort_id.
        """
        with self._get_connection() as conn:
            with conn.cursor() as cur:
                cur.execute(
                    """INSERT INTO cohorts (cohort_id, name, description, project)
                       VALUES (%s, %s, %s, %s)
                       ON CONFLICT (cohort_id) DO UPDATE SET name = EXCLUDED.name;""",
                    (cohort_id, name, description, project),
                )
                if sample_ids:
                    psycopg2.extras.execute_values(
                        cur,
                        "INSERT INTO cohort_members (cohort_id, sample_id) VALUES %s ON CONFLICT DO NOTHING;",
                        [(cohort_id, sid) for sid in sample_ids],
                    )
        logger.info(f"Created cohort '{cohort_id}' with {len(sample_ids)} samples")
        return cohort_id

    def record_qc_metric(
        self,
        sample_id: str,
        pipeline_step: str,
        metric_name: str,
        metric_value: float,
        metric_unit: str = "",
        threshold_min: Optional[float] = None,
        threshold_max: Optional[float] = None,
    ) -> str:
        """Record a QC metric for a sample pipeline step."""
        status = "PASS"
        if threshold_min is not None and metric_value < threshold_min:
            status = "FAIL"
        elif threshold_max is not None and metric_value > threshold_max:
            status = "FAIL"

        sql = """
            INSERT INTO qc_metrics (
                sample_id, pipeline_step, metric_name, metric_value,
                metric_unit, threshold_min, threshold_max, status
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
            RETURNING id;
        """
        with self._get_connection() as conn:
            with conn.cursor() as cur:
                cur.execute(sql, (
                    sample_id, pipeline_step, metric_name, metric_value,
                    metric_unit, threshold_min, threshold_max, status,
                ))
                record_id = str(cur.fetchone()[0])
        return record_id
