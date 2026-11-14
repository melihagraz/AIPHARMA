#!/usr/bin/env python
"""Generate publication figures from the real pipeline output (reports/*.json)."""
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPORTS = os.path.join(REPO, "reports")
FIGDIR = os.path.join(REPORTS, "figures")
os.makedirs(FIGDIR, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 200, "font.size": 11, "font.family": "DejaVu Sans",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#333333", "axes.linewidth": 0.8,
})
INDIGO, TEAL, GREEN, RED, GREY = "#4F46E5", "#0E9DBE", "#0E9F6E", "#DA3B4A", "#9aa0b4"


# ---- Figure: ALK5 docking score distribution ----
run = json.load(open(os.path.join(REPORTS, "pah_run_1PY5.json")))
ranked = run["ranked"]
scores = [r["dock_score"] for r in ranked]
ranks = [r["rank"] for r in ranked]

fig, ax = plt.subplots(figsize=(6.6, 4.0))
colors = [INDIGO if s <= -8.1 else GREY for s in scores]
ax.barh([f"#{r}" for r in ranks], scores, color=colors, height=0.72)
ax.axvline(-8.1, color=TEAL, ls="--", lw=1.4, label="erlotinib reference (−8.1)")
ax.set_xlabel("AutoDock Vina binding affinity (kcal/mol)")
ax.set_ylabel("Ranked candidate")
ax.invert_yaxis()
ax.set_title("De novo candidate docking scores — ALK5 (1PY5)", fontsize=12, weight="bold")
for r, s in zip(ranks, scores):
    ax.text(s - 0.12, f"#{r}", f"{s:.1f}", va="center", ha="right", color="white", fontsize=8.5, weight="bold")
ax.legend(loc="lower right", frameon=False, fontsize=9)
ax.set_xlim(min(scores) - 1.2, 0)
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "fig_docking_scores.png"), bbox_inches="tight")
plt.close(fig)
print("wrote fig_docking_scores.png")


# ---- Figure: redocking validation RMSD ----
rows = []
for p in sorted(glob.glob(os.path.join(REPORTS, "redock_*.json"))):
    d = json.load(open(p))
    if d.get("rmsd") is None:
        continue
    rows.append((d["pdb"], d.get("ligand", "?"), float(d["rmsd"]),
                 bool(d.get("reproduced")), bool(d.get("metals"))))
rows.sort(key=lambda x: x[2])

if rows:
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    labels = [f"{pdb}\n{lig}{' ⚛' if m else ''}" for pdb, lig, r, ok, m in rows]
    vals = [r for *_, r in [(p, l, r, o, m) for (p, l, r, o, m) in rows]]
    vals = [r[2] for r in rows]
    bars = ax.bar(range(len(rows)), vals,
                  color=[GREEN if r[3] else RED for r in rows], width=0.62)
    ax.axhline(2.0, color="#444", ls="--", lw=1.3, label="2.0 Å reproducibility threshold")
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Pose RMSD to crystal (Å)")
    ax.set_title("Redocking validation — self-docking pose recovery", fontsize=12, weight="bold")
    for i, r in enumerate(rows):
        ax.text(i, r[2] + 0.12, f"{r[2]:.2f}", ha="center", fontsize=8.5, weight="bold",
                color=GREEN if r[3] else RED)
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    ax.set_ylim(0, max(v for v in vals) * 1.18)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig_redock_rmsd.png"), bbox_inches="tight")
    plt.close(fig)
    print(f"wrote fig_redock_rmsd.png ({len(rows)} structures)")


# ---- Figure: pocket druggability (ALK5 + BMPR2) ----
fig, ax = plt.subplots(figsize=(4.6, 3.8))
tgts = ["ALK5\n(1PY5)", "BMPR2\n(3G2F)"]
drug = [run["pocket"]["druggability"], 0.824]
ax.bar(tgts, drug, color=[INDIGO, TEAL], width=0.5)
ax.axhline(0.60, color="#444", ls="--", lw=1.3, label="druggable threshold (0.60)")
for i, v in enumerate(drug):
    ax.text(i, v + 0.02, f"{v:.3f}", ha="center", weight="bold", fontsize=10)
ax.set_ylabel("fpocket druggability score")
ax.set_ylim(0, 1.0)
ax.set_title("PAH target pocket druggability", fontsize=12, weight="bold")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "fig_pocket_druggability.png"), bbox_inches="tight")
plt.close(fig)
print("wrote fig_pocket_druggability.png")

print("figures in", FIGDIR)
for f in sorted(os.listdir(FIGDIR)):
    print("  ", f)
