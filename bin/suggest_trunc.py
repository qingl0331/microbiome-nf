#!/usr/bin/env python3
"""
Suggest DADA2 --p-trunc-len-f / --p-trunc-len-r from primer-trimmed paired FASTQs.

Logic (per mate):
  q_pos   = last position before the rolling (window) mean of per-position MEDIAN quality
            first drops below --q
  len_pos = largest length L such that >= --keep-frac of reads are at least L bp long
            (DADA2 discards reads shorter than the truncation length)
  trunc   = min(q_pos, len_pos)
Then enforce merging: trunc_f + trunc_r - amplicon_length >= min_overlap, extending R2 first,
then R1, but never beyond len_pos. If that is impossible the status is FAIL and the
numbers are the best available; review the demux .qzv before using them.

Pure standard library (runs in any container with python3).
"""
import argparse
import gzip
import random
import sys


def profile(files, max_reads, max_files, seed):
    random.seed(seed)
    files = sorted(files)
    if len(files) > max_files:
        files = sorted(random.sample(files, max_files))
    qhist = []      # qhist[pos] = list of 64 counts
    lens = {}
    n = 0
    for fn in files:
        with gzip.open(fn, "rt") as fh:
            i = 0
            while i < max_reads:
                h = fh.readline()
                if not h:
                    break
                fh.readline()
                fh.readline()
                q = fh.readline().rstrip("\n")
                L = len(q)
                lens[L] = lens.get(L, 0) + 1
                while len(qhist) < L:
                    qhist.append([0] * 64)
                for p, ch in enumerate(q):
                    v = ord(ch) - 33
                    if v > 63:
                        v = 63
                    qhist[p][v] += 1
                i += 1
                n += 1
    return qhist, lens, n, len(files)


def median_from_hist(h):
    tot = sum(h)
    if tot == 0:
        return 0
    half = tot / 2.0
    c = 0
    for q, k in enumerate(h):
        c += k
        if c >= half:
            return q
    return 63


def length_cap(lens, n, keep_frac):
    # largest L with fraction(len >= L) >= keep_frac
    at_least = 0
    for L in sorted(lens, reverse=True):
        at_least += lens[L]
        if at_least / n >= keep_frac:
            return L
    return min(lens)


def quality_cap(qhist, q, window):
    med = [median_from_hist(h) for h in qhist]
    for p in range(len(med)):
        lo = max(0, p - window + 1)
        w = med[lo:p + 1]
        if sum(w) / len(w) < q:
            return max(p, 1), med  # p is 0-based index of first failing pos -> truncate before it
    return len(med), med


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--r1", nargs="+", required=True)
    ap.add_argument("--r2", nargs="+", required=True)
    ap.add_argument("--amplicon-length", type=int, required=True,
                    help="expected amplicon length AFTER primer removal (use the upper end of the range)")
    ap.add_argument("--min-overlap", type=int, default=20,
                    help="required overlap; DADA2 needs 12, the rest is safety margin")
    ap.add_argument("--q", type=float, default=25)
    ap.add_argument("--window", type=int, default=10)
    ap.add_argument("--keep-frac", type=float, default=0.95)
    ap.add_argument("--max-reads", type=int, default=20000, help="reads per file")
    ap.add_argument("--max-files", type=int, default=48)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-tsv", required=True)
    a = ap.parse_args()

    res = {}
    for mate, files in (("R1", a.r1), ("R2", a.r2)):
        qhist, lens, n, nf = profile(files, a.max_reads, a.max_files, a.seed)
        if n == 0:
            sys.exit(f"No reads found in {mate} files")
        lcap = length_cap(lens, n, a.keep_frac)
        qcap, med = quality_cap(qhist, a.q, a.window)
        res[mate] = dict(n=n, files=nf, lcap=lcap, qcap=qcap, trunc=min(lcap, qcap), med=med)

    tf, tr = res["R1"]["trunc"], res["R2"]["trunc"]
    need = a.amplicon_length + a.min_overlap
    notes = []
    if tf + tr < need:
        deficit = need - (tf + tr)
        add_r = min(deficit, res["R2"]["lcap"] - tr)
        tr += add_r
        deficit -= add_r
        add_f = min(deficit, res["R1"]["lcap"] - tf)
        tf += add_f
        deficit -= add_f
        notes.append("Quality-based lengths were too short to merge; extended into lower-quality tail "
                     "(DADA2 --p-max-ee still filters bad reads).")
        if deficit > 0:
            notes.append(f"Cannot reach {a.min_overlap} bp overlap for a {a.amplicon_length} bp amplicon even at "
                         f"full read length (short by {deficit} bp). Check --amplicon_length / read length, or "
                         "consider forward-read-only analysis.")
    overlap = tf + tr - a.amplicon_length
    status = "OK" if overlap >= a.min_overlap else "FAIL"

    with open(a.out_tsv, "w") as out:
        out.write("trunc_f\ttrunc_r\texpected_overlap\tstatus\n")
        out.write(f"{tf}\t{tr}\t{overlap}\t{status}\n")

    print("DADA2 truncation suggestion")
    print("===========================")
    for mate in ("R1", "R2"):
        r = res[mate]
        print(f"{mate}: {r['n']} reads from {r['files']} files | quality cap (rolling median >= Q{a.q:g}) = "
              f"{r['qcap']} | length cap ({a.keep_frac:.0%} of reads) = {r['lcap']}")
        meds = r["med"]
        marks = ", ".join(f"{p}:{meds[p-1]}" for p in range(25, len(meds) + 1, 25))
        print(f"    median Q by position -> {marks}")
    print(f"\nSuggested: --trunc_f {tf} --trunc_r {tr}  (expected overlap {overlap} bp on a "
          f"{a.amplicon_length} bp amplicon)  STATUS={status}")
    for n_ in notes:
        print("NOTE: " + n_)


if __name__ == "__main__":
    main()
