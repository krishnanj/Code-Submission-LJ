"""
Part 1 figure — minimal two-panel design.
Panel A: vertical pipeline, outcomes exit right.
Panel B: feature-space scatter, four case types.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

np.random.seed(42)
rng = np.random.default_rng(7)

C = {
    "a": "#27ae60",
    "b": "#7f8c8d",
    "c": "#e67e22",
    "d": "#c0392b",
    "step": "#2c3e50",
}

fig = plt.figure(figsize=(12, 5.5), facecolor="white")
ax1 = fig.add_axes([0.02, 0.06, 0.43, 0.88])
ax2 = fig.add_axes([0.55, 0.12, 0.42, 0.80])


def pill(ax, cx, cy, w, h, label, fc, fs=9.5):
    p = mpatches.FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.10",
        facecolor=fc, edgecolor="none", zorder=3,
    )
    ax.add_patch(p)
    ax.text(cx, cy, label, ha="center", va="center",
            fontsize=fs, color="white", fontweight="bold", zorder=4)


def varr(ax, x, ya, yb, col="#777"):
    ax.annotate("", xy=(x, yb), xytext=(x, ya),
                arrowprops=dict(arrowstyle="-|>", color=col, lw=1.5), zorder=2)


def harr(ax, xa, xb, y, col="#777"):
    ax.annotate("", xy=(xb, y), xytext=(xa, y),
                arrowprops=dict(arrowstyle="-|>", color=col, lw=1.5), zorder=2)


# ─── Panel A ──────────────────────────────────────────────────────────────────
ax1.set_xlim(0, 10); ax1.set_ylim(0, 10); ax1.axis("off")
ax1.text(0.1, 9.85, "A", fontsize=13, fontweight="bold", va="top")

SX = 4.0          # spine center x
BW, BH = 3.6, 0.72

# vertical spine positions
Ys = [8.8, 7.1, 5.4, 3.7, 2.0]

pill(ax1, SX, Ys[0], BW, BH, "Input",            C["step"])
varr(ax1, SX, Ys[0] - BH/2, Ys[1] + BH/2)
pill(ax1, SX, Ys[1], BW, BH, "BLAST",            C["step"])
varr(ax1, SX, Ys[1] - BH/2, Ys[2] + BH/2)
pill(ax1, SX, Ys[2], BW, BH, "FM embedding",     C["step"])
varr(ax1, SX, Ys[2] - BH/2, Ys[3] + BH/2)
pill(ax1, SX, Ys[3], BW, BH, "Ensemble",         C["step"])
varr(ax1, SX, Ys[3] - BH/2, Ys[4] + BH/2, col=C["d"])
pill(ax1, SX, Ys[4], BW, BH, "(d) Adversarial",  C["d"])

# outcome branches (right side)
OX = 8.5; OW = 2.8; OH = 0.66

# from BLAST: (b) no match, (a) confident match — same horizontal line, show via label
harr(ax1, SX + BW/2, OX - OW/2, Ys[1], col=C["b"])
pill(ax1, OX, Ys[1], OW, OH, "(b) Nonsense", C["b"], fs=8.5)
ax1.text(SX + BW/2 + 0.12, Ys[1] + 0.32, "no match", fontsize=7.5,
         color=C["b"], style="italic", ha="left")

# Biological exits at BLAST too — show as note on spine
ax1.text(SX + 0.10, (Ys[0] + Ys[1]) / 2,
         "confident match\n→ (a) Biological",
         fontsize=7.5, color=C["a"], style="italic", ha="left", va="center")

# from FM: (c) far in embedding space
harr(ax1, SX + BW/2, OX - OW/2, Ys[2], col=C["c"])
pill(ax1, OX, Ys[2], OW, OH, "(c) AI-generated", C["c"], fs=8.5)
ax1.text(SX + BW/2 + 0.12, Ys[2] + 0.32, "far in latent space", fontsize=7.5,
         color=C["c"], style="italic", ha="left")

# from Ensemble: (a) Biological (all agree)
harr(ax1, SX + BW/2, OX - OW/2, Ys[3], col=C["a"])
pill(ax1, OX, Ys[3], OW, OH, "(a) Biological", C["a"], fs=8.5)
ax1.text(SX + BW/2 + 0.12, Ys[3] + 0.32, "full agreement", fontsize=7.5,
         color=C["a"], style="italic", ha="left")

ax1.text(SX, 0.8,
         "Krishnan et al. 2026: 6-8 constrained edits can evade\n"
         "DNABERT-2 and Nucleotide Transformer individually",
         ha="center", fontsize=7, color="#aaa", style="italic")


# ─── Panel B ──────────────────────────────────────────────────────────────────
ax2.set_xlim(0, 10); ax2.set_ylim(0, 10)
ax2.set_xlabel("BLAST alignment score", fontsize=10)
ax2.set_ylabel("FM embedding distance", fontsize=10)
ax2.tick_params(labelsize=8)
ax2.set_facecolor("white")

ax2.text(-1.4, 11.2, "B", fontsize=13, fontweight="bold",
         transform=ax2.transData, clip_on=False)
ax2.text(-0.6, 11.2, "Feature space", fontsize=11, color="#333",
         transform=ax2.transData, clip_on=False)

# two simple boundary lines → three regions
VL = 4.2   # BLAST threshold
HL = 5.8   # FM distance threshold

ax2.axvline(VL, color="#bbb", lw=1.3, ls="--", zorder=1)
ax2.axhline(HL, color="#bbb", lw=1.3, ls="--", zorder=1)

# zone shading — simple rectangles
ax2.fill_between([0, VL],    HL, 10, color=C["b"], alpha=0.08)
ax2.fill_between([VL, 10],   HL, 10, color=C["c"], alpha=0.08)
ax2.fill_between([VL, 10], 0, HL,    color=C["a"], alpha=0.08)

# zone labels
ax2.text(2.0, 8.5, "Nonsense", fontsize=9, ha="center",
         color=C["b"], fontstyle="italic")
ax2.text(7.5, 8.5, "AI-generated", fontsize=9, ha="center",
         color=C["c"], fontstyle="italic")
ax2.text(7.5, 2.5, "Biological", fontsize=9, ha="center",
         color=C["a"], fontstyle="italic")

# scatter — 4 case types
bx = rng.uniform(VL + 0.5, 9.5, 22)
by = rng.uniform(0.4, HL - 0.5, 22)
ax2.scatter(bx, by, c=C["a"], s=32, marker="o", alpha=0.85, zorder=5,
            label="(a) Biological")

nx = rng.uniform(0.3, VL - 0.4, 16)
ny = rng.uniform(HL + 0.5, 9.5, 16)
ax2.scatter(nx, ny, c=C["b"], s=32, marker="s", alpha=0.85, zorder=5,
            label="(b) Nonsense")

cx = rng.uniform(VL + 0.5, 9.2, 16)
cy = rng.uniform(HL + 0.5, 9.5, 16)
ax2.scatter(cx, cy, c=C["c"], s=32, marker="^", alpha=0.85, zorder=5,
            label="(c) AI-generated")

# (d) adversarial: starts biological, perturbed close to or across FM boundary
ax2.scatter([7.5], [2.8], c=C["a"], s=32, marker="o", alpha=0.25, zorder=4)
ax2.annotate("", xy=(7.0, HL + 0.4), xytext=(7.5, 3.0),
             arrowprops=dict(arrowstyle="-|>", color="#e67e22", lw=2.0), zorder=6)
ax2.scatter([7.0], [HL + 0.5], c=C["d"], s=140, marker="*",
            edgecolors="#7f0000", linewidths=0.7, zorder=7,
            label="(d) Adversarial")
ax2.text(7.9, HL + 0.4, "6-8 edits", fontsize=8, color="#e67e22",
         fontweight="bold", va="center",
         bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="#e67e22", alpha=0.9))

ax2.legend(fontsize=8.5, loc="lower left", framealpha=0.95,
           markerscale=1.4, handletextpad=0.4, borderpad=0.6, edgecolor="#ccc")

plt.savefig("outputs/figures/part1_figure.png", dpi=150,
            bbox_inches="tight", facecolor="white")
plt.close()
print("Saved.")
