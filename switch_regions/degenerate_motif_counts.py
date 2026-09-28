#!/usr/bin/env python3
# Usage:
#   python3 degenerate_motif_counts.py ighc_upstream5k.fa \
#       -m GAGCW -m GGGBT > degenerate_counts.tsv
#
#   # Also include reverse-complement matches:
#   python3 degenerate_motif_counts.py ighc_upstream5k.fa -m GAGCW -m GGGBT --rc
#
import sys, re, argparse

IUPAC = {
    "A":"A","C":"C","G":"G","T":"T","U":"T","R":"[AG]","Y":"[CT]",
    "S":"[GC]","W":"[AT]","K":"[GT]","M":"[AC]","B":"[CGT]",
    "D":"[AGT]","H":"[ACT]","V":"[ACG]","N":"[ACGT]"
}

def iupac_to_regex(motif: str) -> str:
    out=[]
    for ch in motif.upper():
        out.append(IUPAC.get(ch, re.escape(ch)))
    return "".join(out)

def rc(seq: str) -> str:
    tbl = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")
    return seq.translate(tbl)[::-1].upper().replace("U","T")

def read_fasta(path):
    h=None; b=[]
    with open(path) as f:
        for ln in f:
            ln=ln.strip()
            if not ln: continue
            if ln.startswith(">"):
                if h is not None:
                    yield h, "".join(b).upper().replace("U","T")
                h=ln[1:].strip(); b=[]
            else:
                b.append(ln)
    if h is not None:
        yield h, "".join(b).upper().replace("U","T")

def count_overlapping_regex(seq: str, pat: re.Pattern) -> int:
    # use lookahead to allow overlaps
    return sum(1 for _ in pat.finditer(seq))

def main():
    ap = argparse.ArgumentParser(description="Count degenerate (IUPAC) motifs per FASTA record.")
    ap.add_argument("fasta", help="FASTA file (multi-record ok)")
    ap.add_argument("-m","--motif", action="append", required=True,
                    help="Motif with IUPAC codes (repeat flag for multiple).")
    ap.add_argument("--rc", action="store_true",
                    help="Also count matches on the reverse complement (adds to total).")
    args = ap.parse_args()

    # Compile regexes once (with lookahead)
    regexes = {}
    for m in args.motif:
        rx = iupac_to_regex(m)
        regexes[m] = re.compile(f"(?=({rx}))")

    print("name\tlen_bp\t" +
          "\t".join([f"{m}_count\t{m}_per_kb" for m in args.motif]) +
          "\tGC_fraction")

    for name, seq in read_fasta(args.fasta):
        L = len(seq)
        gc = (seq.count("G")+seq.count("C"))/L if L else 0.0
        fields = [name, str(L)]
        for m, pat in regexes.items():
            c_fwd = count_overlapping_regex(seq, pat)
            if args.rc:
                c_rc = count_overlapping_regex(rc(seq), pat)
            else:
                c_rc = 0
            c = c_fwd + c_rc
            per_kb = c / (L/1000.0) if L else 0.0
            fields += [str(c), f"{per_kb:.2f}"]
        fields += [f"{gc:.3f}"]
        print("\t".join(fields))

if __name__ == "__main__":
    main()
