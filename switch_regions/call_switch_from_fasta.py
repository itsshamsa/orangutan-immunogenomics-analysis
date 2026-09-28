#!/usr/bin/env python3
import sys, os, argparse
from typing import Optional, List
import pandas as pd
import matplotlib.pyplot as plt

# Default motifs typical of IgH switch regions
MOTIFS_DEFAULT = ["GAGCT", "GGGGT"]

def parse_header(h: str):
    """
    Accepts FASTA headers like:
      >IGHG1
      >IGHG1|strand=-
      >IGHG1|chrom=NC_072388.2|start=104863139|end=104868139|strand=-
      >NC_072388.2:104863139-104868139(-)|IGHA1      (also tolerated)
    Returns: (name, chrom, start, end, strand)
    - chrom/start/end are optional and used only to print genomic coords in BED
    - strand defaults to '+'
    """
    name = h
    chrom = None; start = None; end = None; strand = "+"

    # pipe-delimited key=val format
    if "|" in h:
        left, *rest = h.split("|")
        name = left
        kv = {}
        for part in rest:
            if "=" in part:
                k, v = part.split("=", 1)
                kv[k.strip()] = v.strip()
        chrom = kv.get("chrom", chrom)
        if kv.get("start"): 
            try: start = int(kv["start"])
            except: start = None
        if kv.get("end"):   
            try: end = int(kv["end"])
            except: end = None
        strand = kv.get("strand", strand)

    # tolerate UCSC-like chunk "(+)" or "(-)" in the header
    if h.endswith("(+)"): strand = "+"
    if h.endswith("(-)"): strand = "-"

    return name, chrom, start, end, strand

def read_fasta(path: str):
    header = None
    buf = []
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(buf).upper().replace("U","T")
                header = line[1:].strip()
                buf = []
            else:
                buf.append(line)
    if header is not None:
        yield header, "".join(buf).upper().replace("U","T")

def count_overlapping(seq: str, motif: str) -> int:
    """Count overlapping occurrences of motif in seq."""
    i = 0; m = len(motif); cnt = 0
    while True:
        j = seq.find(motif, i)
        if j == -1:
            break
        cnt += 1
        i = j + 1
    return cnt

def sliding_windows(seq: str, win: int, step: int):
    i = 0; n = len(seq)
    while i + win <= n:
        yield i, i + win, seq[i:i+win]
        i += step

def gc_fraction(seq: str) -> float:
    g = seq.count("G"); c = seq.count("C")
    return (g + c) / max(1, len(seq))

def g_skew(seq: str) -> float:
    g = seq.count("G"); c = seq.count("C")
    return (g - c) / max(1, (g + c))

def scan_sequence(name: str, seq: str,
                  chrom: Optional[str], start0: Optional[int], strand: str,
                  win: int, step: int, motifs: List[str],
                  min_gc: float, min_motif_perkb: float, min_len: int):
    """Return (DataFrame of windows, list of (genome_start, genome_end) calls)."""
    import pandas as pd
    rows = []
    for s, e, sub in sliding_windows(seq, win, step):
        gc = gc_fraction(sub)
        gsk = g_skew(sub)
        motif_hits = sum(count_overlapping(sub, m) for m in motifs)
        perkb = motif_hits / (win / 1000.0)
        if chrom is not None and start0 is not None:
            gs = start0 + s
            ge = start0 + e
            out_chrom = chrom
        else:
            # If genome coordinates not provided, write local coords & use sequence name as "chrom"
            gs, ge = s, e
            out_chrom = chrom if chrom else name
        rows.append({
            "seq": name,
            "chrom": out_chrom,
            "start": gs,
            "end": ge,
            "local_start": s,
            "local_end": e,
            "strand": strand,
            "gc_frac": gc,
            "g_skew": gsk,
            "motif_hits": motif_hits,
            "motif_per_kb": perkb
        })
    df = pd.DataFrame(rows)

    # Merge passing windows
    calls = []
    if not df.empty:
        passing = df[(df["gc_frac"] >= min_gc) & (df["motif_per_kb"] >= min_motif_perkb)]
        if not passing.empty:
            passing = passing.sort_values(["local_start", "local_end"])
            cur_s = cur_e = cur_gs = cur_ge = None
            for _, r in passing.iterrows():
                s = int(r["local_start"]); e = int(r["local_end"])
                gs = int(r["start"]);       ge = int(r["end"])
                if cur_s is None:
                    cur_s, cur_e, cur_gs, cur_ge = s, e, gs, ge
                else:
                    if s <= cur_e:  # contiguous/overlapping
                        cur_e = max(cur_e, e)
                        cur_ge = max(cur_ge, ge)
                    else:
                        if (cur_e - cur_s) >= min_len:
                            calls.append((cur_gs, cur_ge))
                        cur_s, cur_e, cur_gs, cur_ge = s, e, gs, ge
            if cur_s is not None and (cur_e - cur_s) >= min_len:
                calls.append((cur_gs, cur_ge))
    return df, calls

def main():
    ap = argparse.ArgumentParser(description="Call switch-like regions from upstream FASTA windows (IgH).")
    ap.add_argument("fasta", help="Input FASTA of upstream windows (use -s with bedtools to preserve strand/orientation)")
    # Defaults below tuned for 5 kb windows. Adjust via CLI if needed.
    ap.add_argument("--win", type=int, default=200, help="Sliding window size (bp)")
    ap.add_argument("--step", type=int, default=50, help="Step size (bp)")
    ap.add_argument("--min-gc", type=float, default=0.55, help="Minimum GC fraction for a window")
    ap.add_argument("--min-motif-perkb", type=float, default=30.0, help="Minimum combined motif density per kb")
    ap.add_argument("--min-len", type=int, default=600, help="Minimum merged length (bp) for a call (use 800–1200 for longer windows)")
    ap.add_argument("--motif", action="append", help="Motif to count (can repeat). Default: GAGCT, GGGGT")
    ap.add_argument("--no-plots", action="store_true", help="Skip QC plots")
    args = ap.parse_args()

    motifs = args.motif if args.motif else MOTIFS_DEFAULT
    base = os.path.splitext(args.fasta)[0]
    bed_out = base + ".switch_calls.bed"
    csv_out = base + ".windows.csv"

    all_rows = []
    bed_lines = []

    for header, seq in read_fasta(args.fasta):
        name, chrom, start, end, strand = parse_header(header)
        df, calls = scan_sequence(
            name, seq, chrom, start, strand,
            args.win, args.step, motifs,
            args.min_gc, args.min_motif_perkb, args.min_len
        )
        all_rows.append(df)
        for gs, ge in calls:
            bed_chrom = chrom if chrom else name
            bed_lines.append(f"{bed_chrom}\t{gs}\t{ge}\t{name}_S_like\t0\t{strand}")

        if not args.no_plots and not df.empty:
            # Simple QC plot: motifs/kb and GC fraction across the window
            x = (df["local_start"] + df["local_end"]) / 2.0
            fig, ax1 = plt.subplots(figsize=(10,3))
            ax1.plot(x, df["motif_per_kb"], label="motifs/kb")
            ax1.set_xlabel("position (bp)")
            ax1.set_ylabel("motifs/kb")
            ax2 = ax1.twinx()
            ax2.plot(x, df["gc_frac"], color="green", label="GC")
            ax2.set_ylabel("GC fraction")
            plt.title(f"{name}")
            fig.tight_layout()
            plt.savefig(f"{base}_{name}_qc.png", dpi=150)
            plt.close(fig)

    if all_rows:
        import pandas as pd
        pd.concat(all_rows, ignore_index=True).to_csv(csv_out, index=False)

    with open(bed_out, "w") as f:
        for line in bed_lines:
            f.write(line + "\n")

    print(f"Wrote: {bed_out}")
    print(f"Wrote: {csv_out}")

if __name__ == "__main__":
    main()
