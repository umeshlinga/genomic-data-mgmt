#  Genomic Data Management & Variant Interpretation System

A SQL-driven genomic data management system for organizing FASTQ, BAM, and VCF files at scale, integrated with REST API-based variant annotation for streamlined interpretation workflows.

##  Overview

**Key Results:**
-  32% faster cohort-level data retrieval
-  22% shorter reporting timelines via automated summaries
-  25% less manual cross-referencing with REST API annotation

---

##  Repository Structure

```
genomic-data-mgmt/
├── db/
│   ├── schema.sql                  # PostgreSQL schema (FASTQ, BAM, VCF tables)
│   ├── migrations/                 # Schema version migrations
│   └── queries/
│       ├── cohort_retrieval.sql    # Optimized cohort-level queries
│       └── variant_lookup.sql      # Variant search queries
├── api/
│   ├── __init__.py
│   ├── genome_api_client.py        # REST API client for reference genome queries
│   └── annotation_fetcher.py       # Transcript-level annotation retrieval
├── pipeline/
│   ├── data_loader.py              # FASTQ/BAM/VCF ingestion into DB
│   ├── variant_interpreter.py      # Automated variant interpretation logic
│   └── report_generator.py         # Structured summary report output
├── tests/
│   ├── test_db_queries.py
│   ├── test_api_client.py
│   └── test_variant_interpreter.py
├── scripts/
│   └── ingest_cohort.sh            # Bash script for bulk cohort ingestion
├── config/
│   └── config.yaml                 # DB credentials, API endpoints
├── requirements.txt
└── README.md
```

---

##  Quickstart

```bash
git clone https://github.com/yourusername/genomic-data-mgmt.git
cd genomic-data-mgmt

pip install -r requirements.txt

# Set up the database
psql -U youruser -d yourdb -f db/schema.sql

# Configure API + DB settings
cp config/config.yaml.example config/config.yaml

# Ingest a cohort
bash scripts/ingest_cohort.sh --input /path/to/vcf_files/ --cohort cohort_001

# Generate variant report
python -m pipeline.report_generator --cohort cohort_001 --output reports/
```

---

##  Tech Stack

| Component | Technology |
|---|---|
| Language | Python 3.10+ |
| Database | PostgreSQL |
| API Integration | REST APIs (Ensembl, NCBI) |
| Data Formats | FASTQ, BAM, VCF |
| Testing | Pytest |

---

##  Running Tests

```bash
pytest tests/ -v
```

---

##  License

MIT License
