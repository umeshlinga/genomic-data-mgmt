"""
genome_api_client.py
--------------------
REST API Client Module

Queries Ensembl REST API and ClinVar for variant annotation,
gene information, and clinical significance data.
"""

import time
import logging
import requests
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ENSEMBL_REST = "https://rest.ensembl.org"
CLINVAR_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@dataclass
class VariantAnnotation:
    """Annotation result for a single variant."""
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str = ""
    transcript_id: str = ""
    consequence: str = ""
    impact: str = ""
    allele_frequency: float = 0.0
    clinvar_significance: str = ""
    protein_change: str = ""
    source: str = "ensembl"

    def to_dict(self) -> dict:
        return {
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "gene": self.gene,
            "transcript_id": self.transcript_id,
            "consequence": self.consequence,
            "impact": self.impact,
            "allele_frequency": self.allele_frequency,
            "clinvar_significance": self.clinvar_significance,
            "protein_change": self.protein_change,
        }


class GenomeAPIClient:
    """
    REST API client for genomic variant annotation.

    Queries Ensembl VEP for functional annotation and
    ClinVar for clinical significance data.

    Parameters
    ----------
    genome_build : str
        Reference genome build: 'GRCh38' or 'GRCh37' (default: 'GRCh38').
    rate_limit_delay : float
        Seconds between API requests to respect rate limits (default: 0.1).
    timeout : int
        Request timeout in seconds (default: 15).
    max_retries : int
        Maximum retry attempts on failure (default: 3).
    """

    def __init__(
        self,
        genome_build: str = "GRCh38",
        rate_limit_delay: float = 0.1,
        timeout: int = 15,
        max_retries: int = 3,
    ):
        self.genome_build = genome_build
        self.rate_limit_delay = rate_limit_delay
        self.timeout = timeout
        self.max_retries = max_retries
        self._session = requests.Session()
        self._session.headers.update({
            "Content-Type": "application/json",
            "Accept": "application/json",
        })

    def _get(self, url: str, params: dict = None) -> Any:
        """Execute a GET request with retry logic."""
        for attempt in range(self.max_retries):
            try:
                time.sleep(self.rate_limit_delay)
                response = self._session.get(url, params=params, timeout=self.timeout)
                response.raise_for_status()
                return response.json()
            except requests.exceptions.HTTPError as e:
                if response.status_code == 429:
                    wait = 2 ** attempt
                    logger.warning(f"Rate limited — waiting {wait}s")
                    time.sleep(wait)
                elif response.status_code == 400:
                    logger.warning(f"Bad request for URL: {url}")
                    return None
                else:
                    raise
            except requests.exceptions.RequestException as e:
                if attempt < self.max_retries - 1:
                    logger.warning(f"Request failed (attempt {attempt+1}): {e}")
                    time.sleep(2 ** attempt)
                else:
                    logger.error(f"All retries failed for {url}: {e}")
                    return None
        return None

    def _post(self, url: str, payload: dict) -> Any:
        """Execute a POST request with retry logic."""
        for attempt in range(self.max_retries):
            try:
                time.sleep(self.rate_limit_delay)
                response = self._session.post(url, json=payload, timeout=self.timeout)
                response.raise_for_status()
                return response.json()
            except requests.exceptions.RequestException as e:
                if attempt < self.max_retries - 1:
                    time.sleep(2 ** attempt)
                else:
                    logger.error(f"POST failed: {e}")
                    return None
        return None

    def annotate_variant(
        self, chrom: str, pos: int, ref: str, alt: str
    ) -> VariantAnnotation:
        """
        Annotate a single variant using Ensembl VEP REST API.

        Parameters
        ----------
        chrom : str
            Chromosome (e.g. 'chr7' or '7').
        pos : int
            Genomic position.
        ref : str
            Reference allele.
        alt : str
            Alternate allele.

        Returns
        -------
        VariantAnnotation
        """
        clean_chrom = chrom.replace("chr", "")
        annotation = VariantAnnotation(chrom=chrom, pos=pos, ref=ref, alt=alt)

        url = f"{ENSEMBL_REST}/vep/human/region/{clean_chrom}:{pos}-{pos}/{alt}"
        params = {
            "content-type": "application/json",
            "canonical": 1,
            "af": 1,
            "af_gnomad": 1,
            "clinvar": 1,
        }

        data = self._get(url, params)
        if not data or not isinstance(data, list):
            return annotation

        result = data[0]
        annotation.consequence = result.get("most_severe_consequence", "")

        # Extract transcript consequences
        tc = result.get("transcript_consequences", [])
        canonical = next((t for t in tc if t.get("canonical") == 1), tc[0] if tc else None)
        if canonical:
            annotation.gene = canonical.get("gene_symbol", "")
            annotation.transcript_id = canonical.get("transcript_id", "")
            annotation.impact = canonical.get("impact", "")
            annotation.protein_change = canonical.get("amino_acids", "")

        # Extract allele frequency from gnomAD
        for coloc in result.get("colocated_variants", []):
            gnomad = coloc.get("gnomad_af") or coloc.get("af", 0)
            if gnomad:
                annotation.allele_frequency = float(gnomad)
                break

        return annotation

    def annotate_batch(
        self, variants: List[Dict], max_batch: int = 200
    ) -> List[VariantAnnotation]:
        """
        Annotate multiple variants using VEP batch endpoint.

        Parameters
        ----------
        variants : list of dict
            Each dict must have keys: chrom, pos, ref, alt.
        max_batch : int
            Maximum variants per API call (default: 200).

        Returns
        -------
        list of VariantAnnotation
        """
        url = f"{ENSEMBL_REST}/vep/human/region"
        annotations = []

        for i in range(0, len(variants), max_batch):
            batch = variants[i:i + max_batch]
            payload = {
                "variants": [
                    f"{v['chrom'].replace('chr','')}:{v['pos']}-{v['pos']}:{v['ref']}/{v['alt']}"
                    for v in batch
                ],
                "canonical": True,
                "af": True,
                "af_gnomad": True,
            }

            results = self._post(url, payload)
            if not results:
                for v in batch:
                    annotations.append(
                        VariantAnnotation(
                            chrom=v["chrom"], pos=v["pos"],
                            ref=v["ref"], alt=v["alt"]
                        )
                    )
                continue

            result_map = {r.get("input", ""): r for r in results}
            for v in batch:
                key = f"{v['chrom'].replace('chr','')}:{v['pos']}-{v['pos']}:{v['ref']}/{v['alt']}"
                result = result_map.get(key, {})
                ann = VariantAnnotation(
                    chrom=v["chrom"], pos=v["pos"],
                    ref=v["ref"], alt=v["alt"]
                )
                if result:
                    ann.consequence = result.get("most_severe_consequence", "")
                    tc = result.get("transcript_consequences", [])
                    canonical = next((t for t in tc if t.get("canonical") == 1), tc[0] if tc else None)
                    if canonical:
                        ann.gene = canonical.get("gene_symbol", "")
                        ann.transcript_id = canonical.get("transcript_id", "")
                        ann.impact = canonical.get("impact", "")
                        ann.protein_change = canonical.get("amino_acids", "")
                annotations.append(ann)

            logger.info(f"Annotated batch {i//max_batch + 1} — {len(batch)} variants")

        return annotations

    def get_gene_info(self, gene_symbol: str) -> Dict:
        """
        Retrieve gene information from Ensembl.

        Parameters
        ----------
        gene_symbol : str
            HGNC gene symbol (e.g. 'BRCA1').

        Returns
        -------
        dict
            Gene metadata including Ensembl ID, location, and description.
        """
        url = f"{ENSEMBL_REST}/lookup/symbol/homo_sapiens/{gene_symbol}"
        params = {"expand": 1, "content-type": "application/json"}
        data = self._get(url, params)
        if not data:
            return {}
        return {
            "ensembl_id": data.get("id", ""),
            "gene_symbol": gene_symbol,
            "chrom": data.get("seq_region_name", ""),
            "start": data.get("start", 0),
            "end": data.get("end", 0),
            "strand": data.get("strand", 0),
            "description": data.get("description", ""),
            "biotype": data.get("biotype", ""),
        }

    def close(self):
        """Close the HTTP session."""
        self._session.close()
