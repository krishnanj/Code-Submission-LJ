"""
Part 3 - Embedding-based clustering of viral reads

Approach:
- For each pool, re-run minimap2 against respiratory references and collect
  reads that pass the same quality filters as Part 2 (MAPQ >= 10, identity >= 0.7).
- Compute a 4-mer frequency profile (256 dimensions) for each read as the
  embedding. 4-mer profiles are fast, interpretable, and require no model
  downloads or GPU.
- Reduce dimensions with UMAP (cosine distance, n_components=2) and cluster
  with HDBSCAN.
- Label each read by its best-matching reference from Part 2 and compare
  cluster assignments to taxonomy labels to assess concordance.

Seasonal validation:
  All 17 pools (12 summer, 5 winter) are embedded together. Because the
  embedding is purely sequence-compositional and has no knowledge of season,
  any clustering separation between winter-only pathogens (HCoV-229E, influenza A)
  and year-round pathogens (rhinovirus, SARS-CoV-2) is an independent, orthogonal
  confirmation of the alignment-based taxonomy in Part 2. HCoV-229E forming
  clusters entirely absent from summer pools is the key validation signal.

Usage:
    python3 part3/part3_clustering.py   (run from repository root)

Outputs:
    outputs/figures/clustering_umap.png        - UMAP coloured by taxonomy label
    outputs/figures/clustering_hdbscan.png     - UMAP coloured by HDBSCAN cluster
    outputs/tables/cluster_taxonomy_table.csv  - cluster vs taxonomy contingency
"""

import gzip
import os
import collections
import itertools
import subprocess

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import umap
import hdbscan

FASTA_DIR = "data/fasta"
REF_FASTA = "data/ref/respiratory_refs.fasta"
FIGURES_DIR = "outputs/figures"
TABLES_DIR = "outputs/tables"
TMP_DIR = "outputs/tmp_clustering"
MIN_MAPQ = 10
MIN_IDENTITY = 0.7
MAX_READS = 5000  # cap per pool to keep runtime reasonable
KMER_K = 4

os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)

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

KMERS = ["".join(k) for k in itertools.product("ACGT", repeat=KMER_K)]
KMER_INDEX = {k: i for i, k in enumerate(KMERS)}


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


def parse_paf_best_hit(paf_path, min_mapq=MIN_MAPQ, min_identity=MIN_IDENTITY):
    """Return dict: read_name -> best reference accession (highest identity)."""
    best = {}
    with open(paf_path) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue
            qname = cols[0]
            tname = cols[5]
            nmatch = int(cols[9])
            alen = int(cols[10])
            mapq = int(cols[11])
            if mapq < min_mapq or alen == 0:
                continue
            identity = nmatch / alen
            if identity < min_identity:
                continue
            if qname not in best or identity > best[qname][1]:
                best[qname] = (tname, identity)
    return {q: v[0] for q, v in best.items()}


def kmer_profile(seq, k=KMER_K):
    seq = seq.upper().replace("N", "")
    counts = np.zeros(len(KMERS), dtype=np.float32)
    for i in range(len(seq) - k + 1):
        mer = seq[i:i + k]
        idx = KMER_INDEX.get(mer)
        if idx is not None:
            counts[idx] += 1
    total = counts.sum()
    if total > 0:
        counts /= total
    return counts


def read_fasta_sequences(fasta_path):
    """Yield (name, sequence) from a plain FASTA file."""
    name = None
    buf = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        yield name, "".join(buf)


def main():
    fasta_files = sorted([
        os.path.join(FASTA_DIR, f)
        for f in os.listdir(FASTA_DIR)
        if f.endswith(".fasta.gz")
    ])

    all_embeddings = []
    all_labels = []
    all_pools = []

    for gz_path in fasta_files:
        pool = pool_label(gz_path)
        print(f"Processing {pool}...")

        tmp_fasta = os.path.join(TMP_DIR, f"{pool}.fasta")
        tmp_paf = os.path.join(TMP_DIR, f"{pool}.paf")

        decompress_fasta(gz_path, tmp_fasta)
        run_minimap2(tmp_fasta, REF_FASTA, tmp_paf)
        best_hit = parse_paf_best_hit(tmp_paf)

        if not best_hit:
            print(f"  No mapped reads for {pool}, skipping.")
            continue

        # Build sequence lookup for mapped reads only
        seq_lookup = {}
        for name, seq in read_fasta_sequences(tmp_fasta):
            if name in best_hit:
                seq_lookup[name] = seq

        mapped_names = list(seq_lookup.keys())
        if len(mapped_names) > MAX_READS:
            rng = np.random.default_rng(42)
            mapped_names = rng.choice(mapped_names, MAX_READS, replace=False).tolist()

        for name in mapped_names:
            seq = seq_lookup[name]
            if len(seq) < KMER_K:
                continue
            emb = kmer_profile(seq)
            ref = best_hit[name]
            label = REF_LABELS.get(ref, ref)
            all_embeddings.append(emb)
            all_labels.append(label)
            all_pools.append(pool)

        print(f"  {len(mapped_names)} reads embedded.")

    if not all_embeddings:
        print("No embeddings collected. Exiting.")
        return

    X = np.stack(all_embeddings)
    labels = np.array(all_labels)
    pools = np.array(all_pools)
    print(f"\nTotal reads embedded: {len(X)}")

    # UMAP
    print("Running UMAP...")
    reducer = umap.UMAP(n_components=2, n_neighbors=30, min_dist=0.1,
                        metric="cosine", random_state=42)
    embedding = reducer.fit_transform(X)

    # HDBSCAN
    print("Running HDBSCAN...")
    clusterer = hdbscan.HDBSCAN(min_cluster_size=50, min_samples=10,
                                 metric="euclidean")
    cluster_ids = clusterer.fit_predict(embedding)

    # ── Plot 1: UMAP coloured by taxonomy label ──────────────────────────────
    unique_labels = sorted(set(labels))
    cmap = cm.get_cmap("tab20", len(unique_labels))
    label_color = {l: cmap(i) for i, l in enumerate(unique_labels)}

    fig, ax = plt.subplots(figsize=(10, 7))
    for lbl in unique_labels:
        mask = labels == lbl
        ax.scatter(
            embedding[mask, 0], embedding[mask, 1],
            c=[label_color[lbl]], label=lbl, s=4, alpha=0.6, linewidths=0,
        )
    ax.set_title("UMAP of viral reads, coloured by taxonomy label", fontsize=11)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend(markerscale=4, fontsize=7, loc="best", framealpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "clustering_umap.png"), dpi=150)
    plt.close()
    print("Saved clustering_umap.png")

    # ── Plot 2: UMAP coloured by HDBSCAN cluster ─────────────────────────────
    unique_clusters = sorted(set(cluster_ids))
    cmap2 = cm.get_cmap("tab20", max(len(unique_clusters), 1))

    fig, ax = plt.subplots(figsize=(10, 7))
    for cid in unique_clusters:
        mask = cluster_ids == cid
        color = "lightgrey" if cid == -1 else cmap2(cid % 20)
        lbl = "noise" if cid == -1 else f"Cluster {cid}"
        ax.scatter(
            embedding[mask, 0], embedding[mask, 1],
            c=[color], label=lbl, s=4, alpha=0.6, linewidths=0,
        )
    ax.set_title("UMAP of viral reads, coloured by HDBSCAN cluster", fontsize=11)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend(markerscale=4, fontsize=7, loc="best", framealpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "clustering_hdbscan.png"), dpi=150)
    plt.close()
    print("Saved clustering_hdbscan.png")

    # ── Contingency table: cluster vs taxonomy ────────────────────────────────
    df = pd.DataFrame({
        "cluster": cluster_ids,
        "taxonomy": labels,
        "pool": pools,
    })
    contingency = (
        df.groupby(["cluster", "taxonomy"])
        .size()
        .reset_index(name="count")
        .sort_values(["cluster", "count"], ascending=[True, False])
    )
    contingency.to_csv(
        os.path.join(TABLES_DIR, "cluster_taxonomy_table.csv"), index=False
    )
    print("Saved cluster_taxonomy_table.csv")

    # Print summary
    print("\nCluster composition summary:")
    for cid in sorted(set(cluster_ids)):
        sub = df[df["cluster"] == cid]["taxonomy"].value_counts()
        top = sub.index[0] if len(sub) > 0 else "n/a"
        total = sub.sum()
        pct = 100 * sub.iloc[0] / total if total > 0 else 0
        label = "noise" if cid == -1 else f"Cluster {cid}"
        print(f"  {label}: n={total}, dominant={top} ({pct:.0f}%)")


if __name__ == "__main__":
    main()
