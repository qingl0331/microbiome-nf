#!/usr/bin/env python3
"""Combine per-step read counts into one wide table (read pairs per sample per step)."""
import argparse, csv, json, os, sys

ap = argparse.ArgumentParser()
ap.add_argument("--fastp", nargs="*", default=[])
ap.add_argument("--counts", nargs="*", default=[])
ap.add_argument("--samples", required=True, help="TSV: sample<TAB>type")
a = ap.parse_args()

types = {}
with open(a.samples) as fh:
    for row in csv.reader(fh, delimiter="\t"):
        if row and row[0] != "sample":
            types[row[0]] = row[1] if len(row) > 1 else "DNA"

tab = {s: {} for s in types}
for fn in a.fastp:
    sid = os.path.basename(fn).rsplit(".fastp.json", 1)[0]
    with open(fn) as fh:
        js = json.load(fh)["summary"]
    # fastp counts reads (both mates); report pairs
    tab.setdefault(sid, {})["fastp_input"] = js["before_filtering"]["total_reads"] // 2
    tab[sid]["after_fastp"] = js["after_filtering"]["total_reads"] // 2
for fn in a.counts:
    with open(fn) as fh:
        for line in fh:
            sid, step, n = line.rstrip("\n").split("\t")
            tab.setdefault(sid, {})[step] = int(n)

steps = ["fastp_input", "after_fastp", "host_removed", "phix_removed", "non_rRNA"]
w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
w.writerow(["sample", "type"] + steps + ["pct_retained"])
for sid in sorted(tab):
    r = tab[sid]
    final = next((r[s] for s in reversed(steps) if s in r), None)
    start = r.get("fastp_input")
    pct = f"{100.0 * final / start:.1f}" if start and final is not None else "NA"
    w.writerow([sid, types.get(sid, "NA")] + [r.get(s, "NA") for s in steps] + [pct])
