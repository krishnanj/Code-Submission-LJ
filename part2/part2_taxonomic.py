"""
Part 2 - Taxonomic classification of Zephyr respiratory virus reads

Approach:
- Map each pool's FASTA file against a curated set of respiratory virus
  reference genomes using minimap2 (splice-unaware, long-read mode).
- Parse PAF output to compute per-reference read counts, genome coverage
  breadth (fraction of reference covered by at least one read), and mean depth.
- Report confidence as mean mapping identity across mapped reads.

Usage:
    python part2_taxonomic.py

Outputs:
    data/results/taxonomy_summary.csv   - per pool per virus summary
    data/results/coverage_heatmap.png   - breadth of coverage heatmap
    data/results/read_counts_heatmap.png - read count heatmap
"""

import subprocess
import gzip
import os
import collections
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

FASTA_DIR = "data/fasta"
REF_FASTA = "data/ref/respiratory_refs.fasta"
RESULTS_DIR = "outputs"
FIGURES_DIR = "outputs/figures"
TABLES_DIR = "outputs/tables"
MIN_MAPQ = 10
MIN_IDENTITY = 0.7

os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)

# ── Reference genome lengths ────────────────────────────────────────────────
def get_ref_lengths(ref_fasta):
    lengths = {}
    current = None
    length = 0
    with open(ref_fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current:
                    lengths[current] = length
                current = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if current:
        lengths[current] = length
    return lengths

# ── Short display names for references ──────────────────────────────────────
REF_LABELS = {
    "NC_001617.1": "Rhinovirus A",
    "NC_001490.1": "Rhinovirus B",
    "NC_009996.1": "Rhinovirus C",
    "NC_038235.1": "RSV A",
    "NC_001781.1": "RSV B",
    "NC_045512.2": "SARS-CoV-2",
    "NC_002645.1": "HCoV-229E",
    "NC_006213.1": "HCoV-OC43",
    "NC_003461.1": "PIV-1",
    "NC_001796.2": "PIV-3",
    "NC_001405.1": "Adenovirus C",
    "NC_004148.2": "Metapneumovirus",
    "NC_004442.1": "Bocaparvovirus",
    "NC_026433.1": "Flu A HA",
    "NC_026435.1": "Flu A PB1",
    "NC_026437.1": "Flu A PA",
    "NC_026438.1": "Flu A PB2",
    "NC_026439.1": "Other",
}

def pool_label(fname):
    base = os.path.basename(fname).replace(".respiratory.fasta.gz", "")
    parts = base.split("-")
    date = parts[0]
    loc = parts[1] if len(parts) > 1 else "?"
    return f"{date[-4:]}-{loc}"

def decompress_fasta(gz_path, out_path):
    with gzip.open(gz_path, "rb") as fin, open(out_path, "wb") as fout:
        fout.write(fin.read())

def run_minimap2(query_fasta, ref_fasta, paf_path):
    cmd = [
        "minimap2",
        "-x", "map-ont",
        "--secondary=no",
        "-c",
        ref_fasta,
        query_fasta,
    ]
    with open(paf_path, "w") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.DEVNULL, check=True)

def parse_paf(paf_path, ref_lengths, min_mapq=MIN_MAPQ, min_identity=MIN_IDENTITY):
    """
    Parse PAF and return per-reference stats.
    PAF columns: qname, qlen, qstart, qend, strand,
                 tname, tlen, tstart, tend, nmatch, alen, mapq, ...
    """
    ref_hits = collections.defaultdict(list)
    with open(paf_path) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue
            tname = cols[5]
            tstart = int(cols[7])
            tend = int(cols[8])
            nmatch = int(cols[9])
            alen = int(cols[10])
            mapq = int(cols[11])
            if mapq < min_mapq or alen == 0:
                continue
            identity = nmatch / alen
            if identity < min_identity:
                continue
            ref_hits[tname].append((tstart, tend, identity))

    results = {}
    for ref, hits in ref_hits.items():
        rlen = ref_lengths.get(ref, 1)
        covered = np.zeros(rlen, dtype=bool)
        for tstart, tend, _ in hits:
            covered[tstart:min(tend, rlen)] = True
        breadth = covered.sum() / rlen
        mean_identity = np.mean([h[2] for h in hits])
        results[ref] = {
            "read_count": len(hits),
            "breadth": breadth,
            "mean_identity": mean_identity,
        }
    return results

def main():
    ref_lengths = get_ref_lengths(REF_FASTA)
    print(f"Loaded {len(ref_lengths)} reference sequences.")

    fasta_files = sorted([
        os.path.join(FASTA_DIR, f)
        for f in os.listdir(FASTA_DIR)
        if f.endswith(".fasta.gz")
    ])
    print(f"Found {len(fasta_files)} pool FASTA files.")

    all_results = {}
    tmp_fasta = os.path.join(TABLES_DIR, "tmp_query.fasta")
    tmp_paf = os.path.join(TABLES_DIR, "tmp_aln.paf")

    for gz_path in fasta_files:
        pool = pool_label(gz_path)
        print(f"Processing pool {pool}...")
        decompress_fasta(gz_path, tmp_fasta)
        run_minimap2(tmp_fasta, REF_FASTA, tmp_paf)
        stats = parse_paf(tmp_paf, ref_lengths)
        all_results[pool] = stats

    os.remove(tmp_fasta)
    os.remove(tmp_paf)

    # ── Write CSV ───────────────────────────────────────────────────────────
    csv_path = os.path.join(TABLES_DIR, "taxonomy_summary.csv")
    refs = sorted(ref_lengths.keys())
    with open(csv_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["pool", "reference", "virus", "read_count",
                         "breadth", "mean_identity"])
        for pool, stats in sorted(all_results.items()):
            for ref in refs:
                s = stats.get(ref, {"read_count": 0, "breadth": 0.0,
                                    "mean_identity": 0.0})
                writer.writerow([
                    pool, ref,
                    REF_LABELS.get(ref, ref),
                    s["read_count"],
                    round(s["breadth"], 4),
                    round(s["mean_identity"], 4),
                ])
    print(f"Saved taxonomy summary to {csv_path}")

    # ── Build dataframes for heatmaps ───────────────────────────────────────
    df = pd.read_csv(csv_path)
    pools = sorted(df["pool"].unique())
    viruses = [REF_LABELS[r] for r in refs if REF_LABELS.get(r, r) != "Other"]
    refs_plot = [r for r in refs if REF_LABELS.get(r, r) != "Other"]

    breadth_mat = pd.DataFrame(index=pools, columns=viruses, dtype=float)
    count_mat = pd.DataFrame(index=pools, columns=viruses, dtype=float)

    for pool in pools:
        for ref in refs_plot:
            virus = REF_LABELS[ref]
            row = df[(df["pool"] == pool) & (df["reference"] == ref)]
            if not row.empty:
                breadth_mat.at[pool, virus] = row["breadth"].values[0]
                count_mat.at[pool, virus] = row["read_count"].values[0]
            else:
                breadth_mat.at[pool, virus] = 0.0
                count_mat.at[pool, virus] = 0.0

    breadth_mat = breadth_mat.astype(float)
    count_mat = count_mat.astype(float)

    # ── Coverage breadth heatmap ─────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(14, 7))
    sns.heatmap(
        breadth_mat,
        ax=ax,
        cmap="YlOrRd",
        vmin=0, vmax=1,
        linewidths=0.4,
        linecolor="grey",
        annot=True,
        fmt=".2f",
        annot_kws={"size": 7},
        cbar_kws={"label": "Genome coverage breadth"},
    )
    ax.set_title("Genome coverage breadth per pool and virus", fontsize=12)
    ax.set_xlabel("Virus")
    ax.set_ylabel("Pool (YYMMDD-location)")
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.yticks(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "coverage_heatmap.png"), dpi=150)
    plt.close()
    print("Saved coverage heatmap.")

    # ── Read count heatmap (log scale) ──────────────────────────────────────
    fig, ax = plt.subplots(figsize=(14, 7))
    log_counts = np.log1p(count_mat.values)
    sns.heatmap(
        pd.DataFrame(log_counts, index=pools, columns=viruses),
        ax=ax,
        cmap="Blues",
        linewidths=0.4,
        linecolor="grey",
        annot=count_mat.values,
        fmt=".0f",
        annot_kws={"size": 7},
        cbar_kws={"label": "log(1 + read count)"},
    )
    ax.set_title("Read counts per pool and virus (colour scale is log)", fontsize=12)
    ax.set_xlabel("Virus")
    ax.set_ylabel("Pool (YYMMDD-location)")
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.yticks(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "read_counts_heatmap.png"), dpi=150)
    plt.close()
    print("Saved read count heatmap.")

    # ── Print top findings ───────────────────────────────────────────────────
    print("\nTop detections (breadth > 0.1):")
    high = df[(df["breadth"] > 0.1) & (df["virus"] != "Other")].sort_values(
        "breadth", ascending=False
    )
    print(high[["pool", "virus", "read_count", "breadth", "mean_identity"]].to_string(index=False))

if __name__ == "__main__":
    main()
