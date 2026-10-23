#!/usr/bin/env python3
"""Characterize the point-lift shift: for reads present in both the mason -oa truth
(reference coords) and the consensus-mapped BAM, compute (lifted cons start) - (truth
start). Tells us whether the 'shifted' fraction is tiny indel-boundary convention
(benign) or a real systematic offset. Run from the validation dir."""
import pysam, bisect, collections, numpy as np

def load(ch):
    s, e, o = [], [], []; t = q = 0
    for L in open(ch):
        L = L.strip()
        if not L:
            continue
        if L.startswith("chain"):
            f = L.split(); t, q = int(f[5]), int(f[10]); continue
        p = L.split(); sz = int(p[0]); s.append(t); e.append(t + sz); o.append(q - t)
        if len(p) == 3:
            t += sz + int(p[1]); q += sz + int(p[2])
    return s, e, o

s, e, o = load("cons2ref.chain")
def lift(x):
    i = bisect.bisect_right(s, x) - 1
    return x + o[i] if i >= 0 and x < e[i] else None

truth = {}
for a in pysam.AlignmentFile("truth.oa.bam", "rb"):
    if a.is_secondary or a.is_supplementary:
        continue
    n = int(a.query_name.rsplit(".", 1)[1])
    truth[n * 2 + (0 if a.is_read1 else 1)] = a.reference_start
print("truth loaded", len(truth), flush=True)

deltas = []
for a in pysam.AlignmentFile("map.cons.bam", "rb"):
    if a.is_secondary or a.is_supplementary or a.is_unmapped:
        continue
    n = int(a.query_name.rsplit(".", 1)[1]); k = n * 2 + (0 if a.is_read1 else 1)
    tp = truth.get(k)
    if tp is None:
        continue
    lp = lift(a.reference_start)
    if lp is None:
        continue
    deltas.append(lp - tp)
    if len(deltas) >= 300000:
        break

d = np.array(deltas)
print("n=%d  exact(0)=%.3f  |d|<=2=%.3f  |d|<=10=%.3f  |d|<=50=%.3f"
      % (len(d), np.mean(d == 0), np.mean(abs(d) <= 2), np.mean(abs(d) <= 10), np.mean(abs(d) <= 50)), flush=True)
print("median|d|=%d mean|d|=%.1f p95|d|=%d max|d|=%d"
      % (np.median(abs(d)), np.mean(abs(d)), np.percentile(abs(d), 95), abs(d).max()), flush=True)
print("delta hist [-5..5]:", dict(sorted(collections.Counter(np.clip(d, -5, 5)).items())), flush=True)
