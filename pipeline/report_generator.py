"""
report_generator.py
-------------------
Variant Interpretation & Report Generator

Generates structured variant interpretation summaries
from annotated VCF data stored in the database.
"""

import json
import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Optional

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Clinical significance priority order
CLINVAR_PRIORITY = {
    "Pathogenic": 1,
    "Likely_pathogenic": 2,
    "Uncertain_significance": 3,
    "Likely_benign": 4,
    "Benign": 5,
}

IMPACT_PRIORITY = {"HIGH": 1, "MODERATE": 2, "LOW": 3, "MODIFIER": 4}


@dataclass
class VariantSummary:
    """Summarized variant for reporting."""
    sample_id: str
    gene: str
    chrom: str
    pos: int
    ref: str
    alt: str
    consequence: str
    impact: str
    allele_frequency: float
    clinvar_significance: str
    protein_change: str
    zygosity: str
    priority_score: int = 0

    def to_dict(self) -> dict:
        return {
            "sample_id": self.sample_id,
            "gene": self.gene,
            "location": f"{self.chrom}:{self.pos}",
            "change": f"{self.ref}>{self.alt}",
            "consequence": self.consequence,
            "impact": self.impact,
            "protein_change": self.protein_change,
            "allele_frequency": self.allele_frequency,
            "zygosity": self.zygosity,
            "clinvar_significance": self.clinvar_significance,
            "priority_score": self.priority_score,
        }


class VariantReportGenerator:
    """
    Generates structured variant interpretation reports.

    Prioritises variants by clinical significance and impact,
    and produces both JSON and HTML summary reports.

    Parameters
    ----------
    output_dir : str
        Directory to write report files.
    """

    def __init__(self, output_dir: str = "reports"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def _compute_priority(self, variant: dict) -> int:
        """Compute a priority score — lower = higher priority."""
        impact_score = IMPACT_PRIORITY.get(variant.get("impact", "MODIFIER"), 4)
        clinvar_score = CLINVAR_PRIORITY.get(
            variant.get("clinvar_significance", ""), 6
        )
        af = variant.get("allele_frequency", 1.0) or 1.0
        af_score = 1 if af < 0.01 else (2 if af < 0.05 else 3)
        return impact_score + clinvar_score + af_score

    def _build_summaries(self, variants: List[Dict]) -> List[VariantSummary]:
        """Convert raw variant dicts to prioritised VariantSummary objects."""
        summaries = []
        for v in variants:
            priority = self._compute_priority(v)
            summaries.append(VariantSummary(
                sample_id=v.get("sample_id", ""),
                gene=v.get("gene", ""),
                chrom=v.get("chrom", ""),
                pos=v.get("pos", 0),
                ref=v.get("ref", ""),
                alt=v.get("alt", ""),
                consequence=v.get("consequence", ""),
                impact=v.get("impact", ""),
                allele_frequency=v.get("allele_frequency", 0.0),
                clinvar_significance=v.get("clinvar_significance", ""),
                protein_change=v.get("protein_change", ""),
                zygosity=v.get("zygosity", ""),
                priority_score=priority,
            ))
        return sorted(summaries, key=lambda x: x.priority_score)

    def generate_json_report(
        self,
        variants: List[Dict],
        sample_id: str,
        metadata: Optional[Dict] = None,
    ) -> str:
        """
        Generate a structured JSON variant report for a sample.

        Parameters
        ----------
        variants : list of dict
            Annotated variant records.
        sample_id : str
            Sample identifier.
        metadata : dict, optional
            Additional metadata to include in the report.

        Returns
        -------
        str
            Path to the JSON report file.
        """
        summaries = self._build_summaries(variants)
        high_impact = [s for s in summaries if s.impact in ("HIGH", "MODERATE")]
        pathogenic = [
            s for s in summaries
            if s.clinvar_significance in ("Pathogenic", "Likely_pathogenic")
        ]

        report = {
            "report_type": "variant_interpretation_summary",
            "generated_at": datetime.now().isoformat(),
            "sample_id": sample_id,
            "metadata": metadata or {},
            "summary": {
                "total_variants": len(variants),
                "high_moderate_impact": len(high_impact),
                "pathogenic_likely_pathogenic": len(pathogenic),
                "unique_genes": len(set(s.gene for s in summaries if s.gene)),
            },
            "top_variants": [s.to_dict() for s in summaries[:50]],
            "pathogenic_variants": [s.to_dict() for s in pathogenic],
        }

        report_path = str(self.output_dir / f"{sample_id}_variant_report.json")
        with open(report_path, "w") as f:
            json.dump(report, f, indent=2)

        logger.info(
            f"JSON report: {report_path} — "
            f"{len(high_impact)} high/moderate, {len(pathogenic)} pathogenic"
        )
        return report_path

    def generate_html_report(
        self,
        variants: List[Dict],
        sample_id: str,
    ) -> str:
        """
        Generate a styled HTML variant interpretation report.

        Parameters
        ----------
        variants : list of dict
        sample_id : str

        Returns
        -------
        str
            Path to the HTML report file.
        """
        summaries = self._build_summaries(variants)[:100]

        impact_colors = {
            "HIGH": "#d32f2f",
            "MODERATE": "#f57c00",
            "LOW": "#388e3c",
            "MODIFIER": "#757575",
        }

        rows = ""
        for s in summaries:
            color = impact_colors.get(s.impact, "#757575")
            af_display = f"{s.allele_frequency:.4f}" if s.allele_frequency else "—"
            rows += f"""
            <tr>
                <td>{s.gene or '—'}</td>
                <td>{s.chrom}:{s.pos}</td>
                <td><code>{s.ref}&gt;{s.alt}</code></td>
                <td>{s.consequence}</td>
                <td style="color:{color};font-weight:bold">{s.impact}</td>
                <td>{s.protein_change or '—'}</td>
                <td>{af_display}</td>
                <td>{s.zygosity}</td>
                <td>{s.clinvar_significance or '—'}</td>
            </tr>"""

        high = sum(1 for s in summaries if s.impact == "HIGH")
        mod = sum(1 for s in summaries if s.impact == "MODERATE")
        path = sum(1 for s in summaries if s.clinvar_significance in ("Pathogenic", "Likely_pathogenic"))

        html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Variant Report — {sample_id}</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 40px; color: #333; }}
  h1 {{ font-size: 22px; }}
  .cards {{ display: flex; gap: 20px; margin: 16px 0 24px; flex-wrap: wrap; }}
  .card {{ background: #f5f5f5; padding: 14px 20px; border-radius: 8px; text-align: center; min-width: 120px; }}
  .card .num {{ font-size: 28px; font-weight: bold; }}
  .card .label {{ font-size: 12px; color: #666; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
  th {{ background: #1565c0; color: white; padding: 9px 12px; text-align: left; }}
  td {{ padding: 7px 12px; border-bottom: 1px solid #eee; }}
  tr:nth-child(even) {{ background: #fafafa; }}
  code {{ background: #f0f0f0; padding: 1px 4px; border-radius: 3px; font-size: 12px; }}
</style>
</head>
<body>
<h1>Variant Interpretation Report</h1>
<p style="color:#666;font-size:13px">Sample: <strong>{sample_id}</strong> &nbsp;|&nbsp; Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
<div class="cards">
  <div class="card"><div class="num">{len(variants)}</div><div class="label">Total variants</div></div>
  <div class="card"><div class="num" style="color:#d32f2f">{high}</div><div class="label">High impact</div></div>
  <div class="card"><div class="num" style="color:#f57c00">{mod}</div><div class="label">Moderate impact</div></div>
  <div class="card"><div class="num" style="color:#7b1fa2">{path}</div><div class="label">Pathogenic</div></div>
</div>
<table>
<thead>
  <tr><th>Gene</th><th>Location</th><th>Change</th><th>Consequence</th>
  <th>Impact</th><th>Protein</th><th>Allele freq</th><th>Zygosity</th><th>ClinVar</th></tr>
</thead>
<tbody>{rows}</tbody>
</table>
</body>
</html>"""

        report_path = str(self.output_dir / f"{sample_id}_variant_report.html")
        with open(report_path, "w") as f:
            f.write(html)

        logger.info(f"HTML report written: {report_path}")
        return report_path

    def generate_cohort_report(
        self,
        cohort_variants: List[Dict],
        cohort_id: str,
    ) -> str:
        """
        Generate a cohort-level gene frequency report.

        Parameters
        ----------
        cohort_variants : list of dict
            Variants across all cohort samples.
        cohort_id : str

        Returns
        -------
        str
            Path to JSON cohort report.
        """
        gene_counts: Dict[str, Dict] = {}
        for v in cohort_variants:
            gene = v.get("gene", "UNKNOWN")
            if gene not in gene_counts:
                gene_counts[gene] = {
                    "gene": gene,
                    "affected_samples": set(),
                    "total_variants": 0,
                    "high_impact": 0,
                    "pathogenic": 0,
                }
            gene_counts[gene]["affected_samples"].add(v.get("sample_id", ""))
            gene_counts[gene]["total_variants"] += 1
            if v.get("impact") == "HIGH":
                gene_counts[gene]["high_impact"] += 1
            if v.get("clinvar_significance") in ("Pathogenic", "Likely_pathogenic"):
                gene_counts[gene]["pathogenic"] += 1

        gene_list = sorted(
            [
                {
                    **{k: v for k, v in g.items() if k != "affected_samples"},
                    "affected_samples": len(g["affected_samples"]),
                }
                for g in gene_counts.values()
            ],
            key=lambda x: (-x["affected_samples"], -x["high_impact"]),
        )

        report = {
            "report_type": "cohort_variant_summary",
            "generated_at": datetime.now().isoformat(),
            "cohort_id": cohort_id,
            "total_variants": len(cohort_variants),
            "unique_genes": len(gene_counts),
            "gene_summary": gene_list[:100],
        }

        report_path = str(self.output_dir / f"{cohort_id}_cohort_report.json")
        with open(report_path, "w") as f:
            json.dump(report, f, indent=2)

        logger.info(f"Cohort report: {report_path} — {len(gene_counts)} genes")
        return report_path
