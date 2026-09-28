#!/usr/bin/env python3
# Usage:
#   python3 compare_species_motifs.py pp_deg_counts.tsv pa_deg_counts.tsv \
#       --out comparison.png
#
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def load_table(path, species):
    df = pd.read_csv(path, sep="\t")
    df["species"] = species
    # simplify gene names (everything after last "|")
    df["gene"] = df["name"].apply(lambda x: x.split("|")[-1].split("::")[0])
    return df

def main():
    ap = argparse.ArgumentParser(description="Compare motif counts between species")
    ap.add_argument("pp", help="pp_deg_counts.tsv (P. pygmaeus)")
    ap.add_argument("pa", help="pa_deg_counts.tsv (P. abelii)")
    ap.add_argument("--out", default="species_comparison.png",
                    help="Output PNG filename")
    args = ap.parse_args()

    df_pp = load_table(args.pp, "P.pygmaeus")
    df_pa = load_table(args.pa, "P.abelii")
    df = pd.concat([df_pp, df_pa])

    # Melt to long form for plotting
    df_long = df.melt(id_vars=["gene","species"], 
                      value_vars=[c for c in df.columns if c.endswith("_per_kb")],
                      var_name="motif", value_name="per_kb")

    # Clean motif names
    df_long["motif"] = df_long["motif"].str.replace("_per_kb","")

    # Sort genes for consistent order
    genes_order = sorted(df_long["gene"].unique())
    motifs = df_long["motif"].unique()

    # Plot
    n_motifs = len(motifs)
    fig, axes = plt.subplots(n_motifs, 1, figsize=(10,4*n_motifs), sharex=True)

    if n_motifs == 1:
        axes = [axes]

    for ax, motif in zip(axes, motifs):
        sub = df_long[df_long["motif"]==motif]
        for sp, color in zip(["P.pygmaeus","P.abelii"], ["steelblue","darkorange"]):
            sub_sp = sub[sub["species"]==sp]
            ax.bar([f"{g}\n" for g in genes_order],
                   sub_sp.set_index("gene").loc[genes_order,"per_kb"],
                   alpha=0.7, label=sp, color=color)
        ax.set_ylabel(f"{motif} / kb")
        ax.set_title(f"{motif} motif density per gene")
        ax.legend()

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(args.out, dpi=150)
    print(f"Wrote {args.out}")

if __name__ == "__main__":
    main()
