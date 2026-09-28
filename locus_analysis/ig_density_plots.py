# ig_density_plots.py
# -------------------
# Creates density summaries and plots from IMGT-style counts and locus lengths.
# Outputs:
#   - CSVs: per-genome rows + species×locus summaries (with 95% CI + pooled density)
#   - Combined PNG plots (bar + scatter)
#   - Separate SVGs per locus (bar + scatter)

import math
from collections import OrderedDict
import pandas as pd
import matplotlib.pyplot as plt

# ============ 1) DATA ============ #
# One row per genome × locus (species, locus, size in bp, total genes)
rows = [
    # Pongo abelii
    ("Pongo abelii", "IGHV", 1088121, 139),
    ("Pongo abelii", "IGHV", 1092016, 140),
    ("Pongo abelii", "IGHV", 1195740, 158),
    ("Pongo abelii", "IGLV",   815129, 104),
    ("Pongo abelii", "IGLV",   810319, 103),
    ("Pongo abelii", "IGLV",   810629, 104),
    ("Pongo abelii", "IGKV",   798671,  77),
    ("Pongo abelii", "IGKV",   806000,  75),
    ("Pongo abelii", "IGKV",   754919,  73),

    # Pongo pygmaeus
    ("Pongo pygmaeus", "IGHV", 1339531, 177),
    ("Pongo pygmaeus", "IGHV", 1270257, 157),
    ("Pongo pygmaeus", "IGLV", 1333666, 113),
    ("Pongo pygmaeus", "IGLV", 1336848, 114),
    ("Pongo pygmaeus", "IGKV",  764914,  73),
    ("Pongo pygmaeus", "IGKV",  806195,  77),

    # Gorilla gorilla gorilla
    ("Gorilla gorilla gorilla", "IGHV",  933797, 135),
    ("Gorilla gorilla gorilla", "IGHV", 1082038, 134),
    ("Gorilla gorilla gorilla", "IGHV",  890210, 120),
    ("Gorilla gorilla gorilla", "IGHV",  893859, 123),
    ("Gorilla gorilla gorilla", "IGLV",  790466,  76),
    ("Gorilla gorilla gorilla", "IGLV",  828625,  86),
    ("Gorilla gorilla gorilla", "IGLV",  781363,  79),
    ("Gorilla gorilla gorilla", "IGLV",  821465,  86),
    ("Gorilla gorilla gorilla", "IGKV",  877825,  43),
    ("Gorilla gorilla gorilla", "IGKV",  671684,  41),
    ("Gorilla gorilla gorilla", "IGKV", 1346898,  50),
    ("Gorilla gorilla gorilla", "IGKV",  866928,  49),

    # Homo sapiens
    ("Homo sapiens", "IGHV",  950000, 170),
    ("Homo sapiens", "IGHV",  960793, 165),
    ("Homo sapiens", "IGHV",  915996, 168),
    ("Homo sapiens", "IGHV",  894048, 172),
    ("Homo sapiens", "IGHV",  983045, 168),
    ("Homo sapiens", "IGLV",  871511,  84),
    ("Homo sapiens", "IGLV",  887529,  95),
    ("Homo sapiens", "IGKV", 1905792,  79),
    ("Homo sapiens", "IGKV", 1166956,  76),
]

df = pd.DataFrame(rows, columns=["Species", "Locus", "Size_bp", "Genes"])

# ============ 2) PER-GENOME DENSITIES ============ #
df["Density_genes_per_Mb"] = df["Genes"] * 1e6 / df["Size_bp"]

# ============ 3) SUMMARY (mean, SD, SE, 95% CI, pooled) ============ #
# Use SciPy t critical if available; else fall back to normal approx 1.96
try:
    from scipy.stats import t as tdist
    def tcrit975(n): return float("nan") if n <= 1 else float(tdist.ppf(0.975, n - 1))
except Exception:
    def tcrit975(n): return float("nan") if n <= 1 else 1.96

grp = df.groupby(["Species", "Locus"])
summary = grp.agg(
    N=("Density_genes_per_Mb", "count"),
    mean_size_bp=("Size_bp", "mean"),
    mean_genes=("Genes", "mean"),
    mean_density=("Density_genes_per_Mb", "mean"),
    sd_density=("Density_genes_per_Mb", "std"),
).reset_index()

summary["SE_density"] = summary["sd_density"] / summary["N"].pow(0.5)
summary["tcrit"] = summary["N"].apply(tcrit975)
summary["CI_low"]  = summary["mean_density"] - summary["tcrit"] * summary["SE_density"]
summary["CI_high"] = summary["mean_density"] + summary["tcrit"] * summary["SE_density"]

# Pooled (length-weighted) density per species×locus
pooled = (df.assign(genes_x1e6=df["Genes"] * 1e6)
            .groupby(["Species", "Locus"])
            .agg(total_bp=("Size_bp", "sum"),
                 total_genes_x1e6=("genes_x1e6", "sum"))
            .assign(pooled_density=lambda x: x["total_genes_x1e6"] / x["total_bp"])
            .reset_index())

summary = summary.merge(pooled[["Species", "Locus", "pooled_density"]],
                        on=["Species", "Locus"], how="left")

# Save CSVs
summary.to_csv("ig_density_summary.csv", index=False)
df.to_csv("ig_per_genome_with_density.csv", index=False)
print("Saved: ig_density_summary.csv, ig_per_genome_with_density.csv")

# ============ 4) COMBINED PLOTS (PNG) ============ #
locus_order = ["IGHV", "IGKV", "IGLV"]
summary["Locus"] = pd.Categorical(summary["Locus"], categories=locus_order, ordered=True)
summary = summary.sort_values(["Locus", "Species"])
species_list = list(OrderedDict.fromkeys(summary["Species"]))
n_species = len(species_list)
bar_width = 0.8 / max(n_species, 1)
x_loci = {l: i for i, l in enumerate(locus_order)}

# Bar: mean density ± SD (all loci/species)
fig, ax = plt.subplots(figsize=(10, 6))
for s_idx, species in enumerate(species_list):
    sub = summary[summary["Species"] == species]
    xs = [x_loci[l] + (s_idx - (n_species-1)/2)*bar_width for l in sub["Locus"]]
    ax.bar(xs, sub["mean_density"], width=bar_width, label=species)
    ax.errorbar(xs, sub["mean_density"], yerr=sub["sd_density"],
                fmt="none", ecolor="black", capsize=4, linewidth=1)
ax.set_xticks(list(x_loci.values()))
ax.set_xticklabels(locus_order)
ax.set_ylabel("Mean density (genes/Mb)")
ax.set_title("Immunoglobulin gene density per species and locus")
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1))
fig.tight_layout()
fig.savefig("ig_density_bar.png", dpi=300)
plt.close(fig)
print("Saved: ig_density_bar.png")
fig.savefig("ig_density_bar.svg", format="svg")
print("Saved: ig_density_bar.svg")


# Scatter: size vs genes per genome
marker_map = {"IGHV": "o", "IGKV": "X", "IGLV": "s"}
fig2, ax2 = plt.subplots(figsize=(9, 6))
for species in df["Species"].unique():
    for locus, marker in marker_map.items():
        sub = df[(df["Species"] == species) & (df["Locus"] == locus)]
        if sub.empty: 
            continue
        ax2.scatter(sub["Size_bp"], sub["Genes"], label=f"{species} – {locus}", marker=marker, s=70)
ax2.set_xlabel("Locus size (bp)")
ax2.set_ylabel("Gene count")
ax2.set_title("Gene count vs. locus size (per genome)")
ax2.legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8)
fig2.tight_layout()
fig2.savefig("ig_size_vs_genes_scatter.png", dpi=300)
plt.close(fig2)
print("Saved: ig_size_vs_genes_scatter.png")
fig2.savefig("ig_size_vs_genes_scatter.svg", format="svg")
print("Saved: ig_size_vs_genes_scatter.svg")


# ============ 5) SEPARATE SVGs PER LOCUS ============ #
# (a) Bar charts per locus
for locus in locus_order:
    sub = summary[summary["Locus"] == locus]
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.bar(sub["Species"], sub["mean_density"])
    ax.errorbar(range(len(sub)), sub["mean_density"], yerr=sub["sd_density"],
                fmt="none", ecolor="black", capsize=4, linewidth=1)
    ax.set_ylabel("Mean density (genes/Mb)")
    ax.set_title(f"Immunoglobulin gene density – {locus}")
    plt.xticks(rotation=30, ha="right")
    fig.tight_layout()
    fname = f"ig_density_{locus}.svg"
    fig.savefig(fname, format="svg")
    plt.close(fig)
    print(f"Saved: {fname}")

# (b) Scatter plots per locus
for locus, marker in marker_map.items():
    sub = df[df["Locus"] == locus]
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for species in sub["Species"].unique():
        ss = sub[sub["Species"] == species]
        ax.scatter(ss["Size_bp"], ss["Genes"], marker=marker, s=70, label=species)
    ax.set_xlabel("Locus size (bp)")
    ax.set_ylabel("Gene count")
    ax.set_title(f"Gene count vs. locus size – {locus}")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8)
    fig.tight_layout()
    fname = f"ig_size_vs_genes_{locus}.svg"
    fig.savefig(fname, format="svg")
    plt.close(fig)
    print(f"Saved: {fname}")
