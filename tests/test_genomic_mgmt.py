"""
test_genomic_mgmt.py
--------------------
Unit tests for the Genomic Data Management System.
"""

import json
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock, call

from pipeline.report_generator import VariantReportGenerator, VariantSummary
from api.genome_api_client import GenomeAPIClient, VariantAnnotation


# ── Fixtures ──────────────────────────────────────────────────────

def make_variant(
    sample_id="SAMPLE_001",
    gene="BRCA1",
    chrom="chr17",
    pos=43092919,
    ref="A",
    alt="T",
    consequence="missense_variant",
    impact="HIGH",
    allele_frequency=0.001,
    clinvar_significance="Pathogenic",
    protein_change="V/L",
    zygosity="HET",
):
    return {
        "sample_id": sample_id,
        "gene": gene,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "consequence": consequence,
        "impact": impact,
        "allele_frequency": allele_frequency,
        "clinvar_significance": clinvar_significance,
        "protein_change": protein_change,
        "zygosity": zygosity,
    }


# ── Report Generator Tests ────────────────────────────────────────

class TestVariantReportGenerator:

    def test_json_report_created(self, tmp_path):
        """JSON report should be created with correct structure."""
        gen = VariantReportGenerator(output_dir=str(tmp_path))
        variants = [
            make_variant(gene="BRCA1", impact="HIGH", clinvar_significance="Pathogenic"),
            make_variant(gene="TP53", impact="MODERATE", clinvar_significance="Uncertain_significance"),
            make_variant(gene="EGFR", impact="LOW", clinvar_significance="Benign"),
        ]
        report_path = gen.generate_json_report(variants, "SAMPLE_001")

        assert Path(report_path).exists()
        with open(report_path) as f:
            report = json.load(f)

        assert report["sample_id"] == "SAMPLE_001"
        assert report["summary"]["total_variants"] == 3
        assert report["summary"]["high_moderate_impact"] == 2
        assert report["summary"]["pathogenic_likely_pathogenic"] == 1
        assert len(report["top_variants"]) == 3

    def test_html_report_created(self, tmp_path):
        """HTML report should contain variant data."""
        gen = VariantReportGenerator(output_dir=str(tmp_path))
        variants = [make_variant(gene="BRCA2", impact="HIGH")]
        report_path = gen.generate_html_report(variants, "SAMPLE_002")

        assert Path(report_path).exists()
        with open(report_path) as f:
            content = f.read()
        assert "SAMPLE_002" in content
        assert "BRCA2" in content
        assert "HIGH" in content

    def test_variant_priority_scoring(self):
        """High impact + pathogenic + rare variant should get lowest priority score."""
        gen = VariantReportGenerator()
        high_rare_path = make_variant(
            impact="HIGH", allele_frequency=0.0001,
            clinvar_significance="Pathogenic"
        )
        low_common_benign = make_variant(
            impact="LOW", allele_frequency=0.3,
            clinvar_significance="Benign"
        )
        score_high = gen._compute_priority(high_rare_path)
        score_low = gen._compute_priority(low_common_benign)
        assert score_high < score_low

    def test_variants_sorted_by_priority(self, tmp_path):
        """Top variants in report should be sorted by priority (highest first)."""
        gen = VariantReportGenerator(output_dir=str(tmp_path))
        variants = [
            make_variant(gene="GENE_LOW", impact="LOW", allele_frequency=0.3,
                        clinvar_significance="Benign"),
            make_variant(gene="GENE_HIGH", impact="HIGH", allele_frequency=0.001,
                        clinvar_significance="Pathogenic"),
            make_variant(gene="GENE_MOD", impact="MODERATE", allele_frequency=0.01,
                        clinvar_significance="Uncertain_significance"),
        ]
        report_path = gen.generate_json_report(variants, "SAMPLE_003")
        with open(report_path) as f:
            report = json.load(f)

        top = report["top_variants"]
        assert top[0]["gene"] == "GENE_HIGH"
        assert top[-1]["gene"] == "GENE_LOW"

    def test_cohort_report(self, tmp_path):
        """Cohort report should aggregate variants by gene correctly."""
        gen = VariantReportGenerator(output_dir=str(tmp_path))
        variants = [
            make_variant(sample_id="S1", gene="TP53", impact="HIGH"),
            make_variant(sample_id="S2", gene="TP53", impact="HIGH"),
            make_variant(sample_id="S1", gene="BRCA1", impact="MODERATE"),
            make_variant(sample_id="S3", gene="EGFR", impact="LOW"),
        ]
        report_path = gen.generate_cohort_report(variants, "COHORT_001")

        with open(report_path) as f:
            report = json.load(f)

        assert report["cohort_id"] == "COHORT_001"
        assert report["total_variants"] == 4
        assert report["unique_genes"] == 3

        tp53 = next(g for g in report["gene_summary"] if g["gene"] == "TP53")
        assert tp53["affected_samples"] == 2
        assert tp53["high_impact"] == 2

    def test_empty_variants_report(self, tmp_path):
        """Empty variant list should produce a valid report with zeros."""
        gen = VariantReportGenerator(output_dir=str(tmp_path))
        report_path = gen.generate_json_report([], "EMPTY_SAMPLE")
        with open(report_path) as f:
            report = json.load(f)
        assert report["summary"]["total_variants"] == 0
        assert report["top_variants"] == []

    def test_metadata_included_in_report(self, tmp_path):
        """Metadata dict should be included in JSON report."""
        gen = VariantReportGenerator(output_dir=str(tmp_path))
        meta = {"project": "TCGA", "library": "WES", "analyst": "Umesh"}
        report_path = gen.generate_json_report(
            [make_variant()], "META_SAMPLE", metadata=meta
        )
        with open(report_path) as f:
            report = json.load(f)
        assert report["metadata"]["project"] == "TCGA"
        assert report["metadata"]["analyst"] == "Umesh"


# ── API Client Tests ──────────────────────────────────────────────

class TestGenomeAPIClient:

    def test_annotation_object_defaults(self):
        """VariantAnnotation should initialise with correct defaults."""
        ann = VariantAnnotation(chrom="chr1", pos=100, ref="A", alt="T")
        assert ann.gene == ""
        assert ann.impact == ""
        assert ann.allele_frequency == 0.0

    def test_annotation_to_dict(self):
        """VariantAnnotation.to_dict() should include all expected keys."""
        ann = VariantAnnotation(
            chrom="chr17", pos=43092919, ref="A", alt="T",
            gene="BRCA1", consequence="missense_variant", impact="HIGH",
            allele_frequency=0.001, clinvar_significance="Pathogenic",
        )
        d = ann.to_dict()
        assert d["gene"] == "BRCA1"
        assert d["consequence"] == "missense_variant"
        assert d["impact"] == "HIGH"
        assert d["allele_frequency"] == 0.001

    @patch("api.genome_api_client.requests.Session.get")
    def test_annotate_variant_success(self, mock_get):
        """Successful API call should return populated VariantAnnotation."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = [{
            "most_severe_consequence": "missense_variant",
            "transcript_consequences": [{
                "gene_symbol": "BRCA1",
                "transcript_id": "ENST00000357654",
                "impact": "HIGH",
                "canonical": 1,
                "amino_acids": "V/L",
            }],
            "colocated_variants": [{"gnomad_af": 0.0005}],
        }]
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        client = GenomeAPIClient()
        ann = client.annotate_variant("chr17", 43092919, "A", "T")

        assert ann.gene == "BRCA1"
        assert ann.consequence == "missense_variant"
        assert ann.impact == "HIGH"
        assert ann.allele_frequency == 0.0005
        assert ann.transcript_id == "ENST00000357654"

    @patch("api.genome_api_client.requests.Session.get")
    def test_annotate_variant_api_failure(self, mock_get):
        """API failure should return empty VariantAnnotation without raising."""
        mock_get.side_effect = Exception("Connection refused")
        client = GenomeAPIClient(max_retries=1)
        ann = client.annotate_variant("chr1", 100, "A", "T")
        assert ann.gene == ""
        assert ann.consequence == ""

    @patch("api.genome_api_client.requests.Session.get")
    def test_get_gene_info(self, mock_get):
        """Gene info lookup should return correctly structured dict."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "id": "ENSG00000012048",
            "seq_region_name": "17",
            "start": 43044295,
            "end": 43125483,
            "strand": -1,
            "description": "BRCA1 DNA repair associated",
            "biotype": "protein_coding",
        }
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        client = GenomeAPIClient()
        info = client.get_gene_info("BRCA1")

        assert info["ensembl_id"] == "ENSG00000012048"
        assert info["chrom"] == "17"
        assert info["biotype"] == "protein_coding"
