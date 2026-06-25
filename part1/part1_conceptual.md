# Part 1 — Evaluation Design: AI-Generated vs Biological Nucleic Acid Sequences

An incoming sequence falls into one of four categories (Figure 1A). **(a) Biological:** BLAST returns a confident high-identity match to a known species. **(b) Nonsense:** BLAST finds nothing; flagged immediately. **(c) AI-generated:** BLAST returns a partial match, but the sequence is distant from real references in the latent space of a genomic foundation model (FM) such as DNABERT-2 or the Nucleotide Transformer. Cases (a)-(c) are handled by pairing a deterministic method (BLAST) with a probabilistic one (FM embedding comparison).

**(d) Adversarially perturbed:** the hard case. Biologically constrained edits keep both BLAST identity and FM embedding similarity high while meaningfully altering the sequence. Krishnan et al. (2026) showed 6-8 such edits suffice to evade DNABERT-2 and the Nucleotide Transformer individually. Case (d) requires an ensemble pipeline combining BLAST scores, FM distances, codon usage deviation, and motif conservation, where cross-method disagreement is itself a flag (Figure 1B).

The benchmark covers all four categories, with category (d) constructed using the constrained-edit procedure from Krishnan et al. (2026). Metrics: AUROC and precision-recall per category, plus calibration — reliable confidence scores matter as much as rank ordering in biosurveillance.

---

![Figure 1](data/results/part1_figure.png)

*Figure 1. (A) Screening pipeline for the four sequence categories. (B) Feature space showing how cases map onto BLAST score vs FM embedding distance, and why adversarially perturbed sequences require ensemble screening to distinguish from genuine biological sequences.*

---

**Reference**

Krishnan J, Rangarajan AM, Loehr A, Hoelscher-Obermaier J. Adversarial Genomic Sequences Could Evade Biosecurity Screening. Presented at CyberBio 2026 (workshop collocated with IEEE Symposium on Security and Privacy 2026); best paper runner-up; in press, IEEE. Also accepted at ICLR Workshop on Machine Learning for Genomics Explorations (MLGenX) 2026. https://openreview.net/forum?id=j0o93ydeGW
