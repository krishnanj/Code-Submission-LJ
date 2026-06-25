# Zephyr Metagenomic Biosurveillance — Analysis Pipeline

**Jeyashree Krishnan · June 2026**

Analysis pipeline for SecureBio Zephyr metagenomic biosurveillance data: taxonomic classification of 17 pooled nasal-swab sequencing runs (12 summer, 5 winter) and unsupervised embedding-based clustering of viral reads.

---

## Repository structure

```
Code-Submission-LJ/
├── part1/
│   ├── part1_conceptual.md     # Written framework: detecting AI-generated sequences
│   └── part1_figure.py         # Generates outputs/figures/part1_figure.png
│
├── part2/
│   └── part2_taxonomic.py      # minimap2-based taxonomic classification (all 17 pools)
│
├── part3/
│   └── part3_clustering.py     # 4-mer embedding → UMAP → HDBSCAN clustering
│
├── data/
│   ├── ref/
│   │   └── respiratory_refs.fasta   # 17 curated respiratory virus reference genomes
│   └── fasta/                       # FASTA.gz pool files (not committed — see below)
│
└── outputs/
    ├── figures/
    │   ├── part1_figure.png          # Detection framework diagram
    │   ├── coverage_heatmap.png      # Genome breadth across all 17 pools × 17 references
    │   ├── read_counts_heatmap.png   # Read counts (log scale)
    │   ├── clustering_umap.png       # UMAP projection coloured by species
    │   └── clustering_hdbscan.png    # HDBSCAN cluster assignments
    └── tables/
        ├── taxonomy_summary.csv      # Full per-pool per-species results
        └── cluster_taxonomy_table.csv # Cluster ID → dominant taxonomy
```

---

## Dependencies

**Python packages**
```bash
pip install numpy pandas matplotlib seaborn umap-learn hdbscan
```

**minimap2** (must be on PATH)
```bash
brew install minimap2        # macOS
# or: conda install -c bioconda minimap2
```

---

## Data

FASTA files are not committed due to size. Download from the [SecureBio Zephyr public data portal](https://securebio.org/zephyr).

| Cohort | Pools | Date range | Pool ID prefix |
|--------|-------|------------|----------------|
| Summer | 12 | May–June 2026 | `26MMDD-*` |
| Winter (seasonal validation) | 5 | Nov–Dec 2025 | `25MMDD-*` |

The 17 pool IDs analyzed are listed in `outputs/tables/taxonomy_summary.csv`.

Place downloaded `.respiratory.fasta.gz` files in `data/fasta/` before running Parts 2 or 3.

---

## Running the analysis

All scripts are run from the **repository root**.

### Part 1 — detection framework figure
```bash
python3 part1/part1_figure.py
# Output: outputs/figures/part1_figure.png
```

### Part 2 — taxonomic classification
```bash
python3 part2/part2_taxonomic.py
# Outputs: outputs/tables/taxonomy_summary.csv
#          outputs/figures/coverage_heatmap.png
#          outputs/figures/read_counts_heatmap.png
# Runtime: ~30–60 min (minimap2 on ONT reads, 17 pools)
```

### Part 3 — embedding-based clustering
```bash
python3 part3/part3_clustering.py
# Outputs: outputs/figures/clustering_umap.png
#          outputs/figures/clustering_hdbscan.png
#          outputs/tables/cluster_taxonomy_table.csv
# Runtime: ~15–30 min (UMAP + HDBSCAN on ~60k reads)
```

---

## Methods

### Part 2 — Taxonomic classification

| Parameter | Value |
|-----------|-------|
| Aligner | minimap2 `-x map-ont --secondary=no -c` |
| MAPQ filter | ≥ 10 |
| Identity filter | ≥ 0.70 |
| Primary metric | **Genome coverage breadth** (fraction of reference covered by ≥ 1 read) |
| References | 17 genomes: Rhinovirus A/B/C, RSV A/B, SARS-CoV-2, HCoV-229E, HCoV-OC43, PIV-1/3, Adenovirus C, Metapneumovirus, Bocaparvovirus, Influenza A (HA, PB1, PA, PB2 as separate segments) |

Breadth is the primary confidence metric because it distinguishes genuine infections from cross-mapping artifacts more reliably than read depth alone.

### Part 3 — Embedding-based clustering

| Parameter | Value |
|-----------|-------|
| Embedding | 4-mer frequency profiles, 256 dimensions, L1-normalised per read |
| Subsampling | Up to 5,000 reads per pool (seed 42) |
| Dim. reduction | UMAP — cosine distance, `n_neighbors=30`, `min_dist=0.1` |
| Clustering | HDBSCAN — `min_cluster_size=50`, `min_samples=10` |
| Labels | Best-hit accession from Part 2 |

---

## Key findings

- **Rhinovirus A and B** co-circulate year-round across all sampled Boston sites (100% genome breadth in multiple pools both seasons).
- **SARS-CoV-2** detected in both seasons; strongest signal in winter pool 1206-Copl (100% breadth, 102,983 reads).
- **HCoV-229E** and **Influenza A** are exclusively winter detections (Nov–Dec 2025 cohort), consistent with temperate-climate seasonality — the pipeline recovers this pattern without any seasonal tuning.
- Influenza A confirmed by multi-segment concordance: PA, PB1, PB2, and HA all above 70% breadth in Dec 2025 pools.
- Unsupervised 4-mer clustering (490 clusters, 56,833 summer reads) independently validates alignment taxonomy. Rhinovirus A and B resolve into hundreds of subclusters reflecting serotype-level diversity. SARS-CoV-2 and HCoV-OC43 form ~17 and ~25 tight, well-separated clusters respectively.

---

## Reference

Krishnan J, Rangarajan AM, Loehr A, Hoelscher-Obermaier J. Adversarial Genomic Sequences Could Evade Biosecurity Screening. Presented at CyberBio 2026 (workshop collocated with IEEE Symposium on Security and Privacy 2026); best paper runner-up; in press, IEEE. Also accepted at ICLR Workshop on Machine Learning for Genomics Explorations (MLGenX) 2026. https://openreview.net/forum?id=j0o93ydeGW
