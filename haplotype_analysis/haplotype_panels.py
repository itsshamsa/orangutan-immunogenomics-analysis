import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# -------- settings --------
excel_path = "haplotype PA.xlsx"   # change if needed
LOCI = ["IGH", "IGK", "IGL"]       # your columns
SAVE_FMT = "svg"
# --------------------------

def read_and_pair(path):
    df_raw = pd.read_excel(path)
    cols = list(df_raw.columns)

    # map each locus to (H1, H2) = (its column, the very next column)
    pairs = {}
    for locus in LOCI:
        try:
            i = cols.index(locus)
            h1 = cols[i]
            h2 = cols[i + 1]  # the next column (often "Unnamed: X")
            pairs[locus] = (h1, h2)
        except (ValueError, IndexError):
            pass

    missing = [L for L in LOCI if L not in pairs]
    if missing:
        raise SystemExit(f"Could not detect both haplotype columns for: {', '.join(missing)}.\nColumns seen: {cols}")

    # rename to standard names
    ren = {}
    for L in LOCI:
        c1, c2 = pairs[L]
        ren[c1] = f"{L}_H1"
        ren[c2] = f"{L}_H2"
    df = df_raw.rename(columns=ren)

    # drop the first row if it looks like hap1/hap2 labels
    first = df.iloc[0]
    lab_vals = []
    for L in LOCI:
        for side in ("H1", "H2"):
            col = f"{L}_{side}"
            if col in df.columns:
                v = first[col]
                if pd.notna(v):
                    lab_vals.append(str(v).strip().lower())
    if any(v in ("hap1", "hap 1", "hap2", "hap 2") for v in lab_vals):
        df = df.iloc[1:].reset_index(drop=True)

    # keep only standardized six columns
    ordered = []
    for L in LOCI:
        ordered += [f"{L}_H1", f"{L}_H2"]
    return df[ordered]

def tidy_for_locus(df, locus):
    sub = df[[f"{locus}_H1", f"{locus}_H2"]].rename(columns={f"{locus}_H1": "H1", f"{locus}_H2": "H2"})
    sub = sub.dropna(how="all")

    # gene label = part before '*' from whichever side is present
    def label_row(r):
        for v in (r["H1"], r["H2"]):
            if pd.notna(v) and str(v).strip():
                return str(v).split("*")[0]
        return "NA"

    labels = sub.apply(label_row, axis=1)
    counts, uniq = {}, []
    for g in labels:
        counts[g] = counts.get(g, 0) + 1
        uniq.append(g if counts[g] == 1 else f"{g}({counts[g]})")
    sub["gene_label"] = uniq
    return sub

def statuses(h1, h2):
    h1_na = h1.isna() | (h1.astype(str).str.strip() == "")
    h2_na = h2.isna() | (h2.astype(str).str.strip() == "")
    same  = (~h1_na) & (~h2_na) & (h1.astype(str) == h2.astype(str))
    # Any missing (either side) -> red
    h1_status = np.where(h1_na, "missing", "present")
    h2_status = np.where(h2_na, "missing", np.where(same, "same", "different"))
    return h1_status, h2_status

COLORS = {
    "same": "#377eb8",      # blue
    "present": "#377eb8",   # blue (H1 present)
    "different": "orange",  # orange
    "missing": "red",       # red for ANY missing
}

df = read_and_pair(excel_path)

for locus in LOCI:
    sub = tidy_for_locus(df, locus)
    if sub.empty:
        continue

    s1, s2 = statuses(sub["H1"], sub["H2"])
    x = np.arange(len(sub))

    plt.figure(figsize=(16, 3.2))
    # H1 on TOP (y=1), H2 on BOTTOM (y=0)
    for i in range(len(sub)):
        plt.bar(x[i], 0.9, bottom=1, width=0.95,
                color=COLORS[s1[i]], edgecolor="black", linewidth=0.2)  # H2 top
        plt.bar(x[i], 0.9, bottom=0, width=0.95,
                color=COLORS[s2[i]], edgecolor="black", linewidth=0.2)  # H1 bottom

    plt.xlim(-0.5, len(sub)-0.5)
    plt.ylim(-0.1, 2.1)
    plt.yticks([1.45, 0.45], ["Haplotype 2", "Haplotype 1"])
    plt.xticks(x, sub["gene_label"], rotation=75, ha="right", fontsize=7)
    plt.title(locus, fontweight="bold")
    plt.xlabel(f"{locus} gene")

    legend_elems = [
        Patch(fc=COLORS["same"], ec="black", label="Same / Present"),
        Patch(fc=COLORS["different"], ec="black", label="Different"),
        Patch(fc=COLORS["missing"], ec="black", label="Missing"),
    ]
    plt.legend(handles=legend_elems, ncol=3, loc="upper left", frameon=False)

    plt.tight_layout()
    out_path = f"{locus}_haplotype_panel.{SAVE_FMT}"
    plt.savefig(out_path, format=SAVE_FMT, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_path}")
