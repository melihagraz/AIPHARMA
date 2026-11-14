#!/usr/bin/env python
"""Redocking pose overlay (3D) + pipeline diagram figures."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(REPO, "reports", "figures")
os.makedirs(FIGDIR, exist_ok=True)
INDIGO, TEAL, GREEN, MAG = "#4F46E5", "#0E9DBE", "#16a34a", "#db2777"


def read_atoms(path):
    xs = []
    for l in open(path):
        if l.startswith(("ATOM", "HETATM")):
            xs.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(xs)


def bond_list(coords, cutoff=1.95):
    out = []
    n = len(coords)
    for i in range(n):
        for j in range(i + 1, n):
            if np.linalg.norm(coords[i] - coords[j]) < cutoff:
                out.append((i, j))
    return out


# ---- Figure: 3D redocking pose overlay (1PY5) ----
cryst = read_atoms(os.path.join(REPO, "reports", "redock_1PY5_crystal.pdb"))
dock = read_atoms(os.path.join(REPO, "reports", "redock_1PY5_docked.pdb"))

fig = plt.figure(figsize=(6.2, 5.4))
ax = fig.add_subplot(111, projection="3d")
for coords, col, lab in [(cryst, GREEN, "crystal pose (experimental)"),
                         (dock, MAG, "re-docked pose")]:
    for i, j in bond_list(coords):
        ax.plot(*zip(coords[i], coords[j]), color=col, lw=3.2, solid_capstyle="round")
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], color=col, s=26, depthshade=True)
    ax.plot([], [], color=col, lw=3, label=lab)

ax.set_title("Redocking pose overlay — ALK5 (1PY5, ligand PY1)\nRMSD 0.46 Å  ·  pose reproduced",
             fontsize=11.5, weight="bold")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0.03), frameon=False, ncol=2, fontsize=9.5)
ax.set_axis_off()
ax.view_init(elev=18, azim=-60)
try:
    ax.set_box_aspect((1, 1, 1))
except Exception:
    pass
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "fig_redock_overlay.png"), dpi=200, bbox_inches="tight")
plt.close(fig)
print("wrote fig_redock_overlay.png")


# ---- Figure: pipeline / graphical abstract ----
stages = [
    ("Target\nstructure", "PDB / RCSB", "#e9e7fb"),
    ("Pocket\ndetection", "fpocket\ndruggability", "#dbeafe"),
    ("Pharmacophore", "RDKit\n3D features", "#d1fae5"),
    ("De novo\ndesign", "REINVENT", "#fef3c7"),
    ("ADMET\nfilter", "ADMET-AI\nAmes gate", "#fee2e2"),
    ("Docking +\nranking", "AutoDock Vina", "#ede9fe"),
]
fig, ax = plt.subplots(figsize=(11, 2.7))
ax.set_xlim(0, len(stages) * 2.05)
ax.set_ylim(0, 2)
x = 0.15
for i, (title, sub, col) in enumerate(stages):
    box = FancyBboxPatch((x, 0.55), 1.55, 0.95, boxstyle="round,pad=0.04,rounding_size=0.12",
                         linewidth=1.1, edgecolor="#4b4b6b", facecolor=col)
    ax.add_patch(box)
    ax.text(x + 0.775, 1.18, title, ha="center", va="center", fontsize=10.5, weight="bold", color="#1b1840")
    ax.text(x + 0.775, 0.82, sub, ha="center", va="center", fontsize=8.2, color="#3a3560")
    if i < len(stages) - 1:
        ax.add_patch(FancyArrowPatch((x + 1.58, 1.02), (x + 2.02, 1.02),
                     arrowstyle="-|>", mutation_scale=16, color="#4F46E5", lw=1.6))
    x += 2.05
ax.text(len(stages) * 1.025, 0.2, "reproducibility validated by self-docking (RMSD < 2 Å); metalloprotein-aware receptor prep",
        ha="center", fontsize=8.4, style="italic", color="#565080")
ax.set_axis_off()
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "fig_pipeline.png"), dpi=200, bbox_inches="tight")
plt.close(fig)
print("wrote fig_pipeline.png")
