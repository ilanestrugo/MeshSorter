#!/usr/bin/env python3
"""What does the order-derived destination sequence look like? (Section S9)

No simulation is involved.  The script reads results/order_derived_sequence.txt,
the belt label of each of the 98,816 usable orders in chronological order, and
reports:

  1. the share of each belt;
  2. for lags 1 to MAXLAG, the chi-square test of independence between the label
     of an order and the label of the order `lag` places later, with Cramer's V,
     so that a small dependence at lag one is not taken for weak dependence at
     every lag;
  3. what each feeder sees when the sequence is consumed by m feeders in index
     order, which is how the simulator uses it: feeder j reads positions
     j, j+m, j+2m, ... so each feeder's stream is a fixed subsequence of
     N/m labels, repeated cyclically;
  4. the same lag statistics for random permutations of the sequence, which have
     the same counts and no chronology, so that the values above can be read
     against what chance alone produces.

    python3 sequence_dependence.py
"""
import csv, json, math, os, random, sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
SEQ = os.path.join(OUT, "order_derived_sequence.txt")
BELTS, FEEDERS, MAXLAG = 4, 4, 20
PERMS = 200

x = np.array([int(l) for l in open(SEQ) if l.strip() and not l.startswith("#")]) - 1
N = len(x)


def lag_table(seq, lag):
    t = np.zeros((BELTS, BELTS))
    np.add.at(t, (seq[:-lag], seq[lag:]), 1)
    return t


def lag_stat(seq, lag):
    t = lag_table(seq, lag)
    chi2, p, dof, _ = stats.chi2_contingency(t, correction=False)
    v = math.sqrt(chi2 / (t.sum() * (BELTS - 1)))
    return chi2, p, v


res = {}
counts = np.bincount(x, minlength=BELTS)
res["labels"] = int(N)
res["counts"] = counts.tolist()
res["shares"] = (counts / N).tolist()
print(f"{N} orders; belt counts {counts.tolist()}")
print("shares     " + ", ".join(f"{s:.4f}" for s in counts / N))
print(f"expected under a uniform draw: 0.2500 each; largest departure "
      f"{np.max(np.abs(counts / N - 0.25)):.4f}")
chi2, p = stats.chisquare(counts)
print(f"chi-square goodness of fit to equal shares: chi2 = {chi2:.2f}, df = 3, p = {p:.3f}")
res["share_gof"] = dict(chi2=float(chi2), p=float(p))

print(f"\nlag   chi2     p        Cramer's V")
lags = []
for lag in range(1, MAXLAG + 1):
    c, pv, v = lag_stat(x, lag)
    lags.append(dict(lag=lag, chi2=float(c), p=float(pv), V=float(v)))
    print(f"{lag:3d}  {c:7.2f}  {pv:7.4f}  {v:.4f}")
res["lags"] = lags
nsig = sum(1 for r in lags if r["p"] < 0.05)
nbonf = sum(1 for r in lags if r["p"] < 0.05 / MAXLAG)
print(f"\nlags with p < 0.05: {nsig} of {MAXLAG};  with p < 0.05/{MAXLAG}: {nbonf}")
res["lags_sig_05"], res["lags_sig_bonferroni"] = nsig, nbonf

# what chance alone gives for these statistics
rng = random.Random(20261004)
vmax1, vany = [], []
for _ in range(PERMS):
    perm = x.copy()
    np.random.RandomState(rng.randrange(2 ** 31)).shuffle(perm)
    vs = [lag_stat(perm, lag)[2] for lag in range(1, MAXLAG + 1)]
    vmax1.append(vs[0])
    vany.append(max(vs))
v1 = lags[0]["V"]
vmax_obs = max(r["V"] for r in lags)
print(f"\nrandom permutations of the same labels ({PERMS}):")
print(f"  lag-1 V: mean {np.mean(vmax1):.4f}, 95th percentile {np.percentile(vmax1, 95):.4f}; "
      f"observed {v1:.4f}")
print(f"  largest V over lags 1..{MAXLAG}: mean {np.mean(vany):.4f}, "
      f"95th percentile {np.percentile(vany, 95):.4f}; observed {vmax_obs:.4f}")
res["permutation"] = dict(perms=PERMS, lag1_V_mean=float(np.mean(vmax1)),
                          lag1_V_p95=float(np.percentile(vmax1, 95)),
                          maxV_mean=float(np.mean(vany)),
                          maxV_p95=float(np.percentile(vany, 95)),
                          observed_lag1_V=v1, observed_maxV=float(vmax_obs))

# slow drift: pairwise lag tables cannot see it
BLOCKS = 98
size = N // BLOCKS
tab = np.array([np.bincount(x[b * size:(b + 1) * size], minlength=BELTS) for b in range(BLOCKS)],
               dtype=float)
c, pv, dof, _ = stats.chi2_contingency(tab, correction=False)
null = []
for _ in range(PERMS):
    perm = x.copy()
    np.random.RandomState(rng.randrange(2 ** 31)).shuffle(perm)
    t2 = np.array([np.bincount(perm[b * size:(b + 1) * size], minlength=BELTS)
                   for b in range(BLOCKS)], dtype=float)
    null.append(stats.chi2_contingency(t2, correction=False)[0])
print(f"\nslow drift: {BLOCKS} consecutive blocks of {size} orders, belt shares compared across blocks")
print(f"  chi-square = {c:.1f} on {dof} df, p = {pv:.3f}; "
      f"random permutations give mean {np.mean(null):.1f}, 95th percentile {np.percentile(null, 95):.1f}")
res["blocks"] = dict(blocks=BLOCKS, size=int(size), chi2=float(c), dof=int(dof), p=float(pv),
                     perm_mean=float(np.mean(null)), perm_p95=float(np.percentile(null, 95)))

# feeder view
print(f"\nthe sequence as {FEEDERS} feeders consume it (feeder j reads positions j, j+{FEEDERS}, ...)")
print(f"N = {N} = {FEEDERS} x {N // FEEDERS}" + ("" if N % FEEDERS == 0 else f" + {N % FEEDERS}"))
feeder = []
for j in range(FEEDERS):
    sub = x[j::FEEDERS]
    sh = np.bincount(sub, minlength=BELTS) / len(sub)
    feeder.append(dict(feeder=j + 1, length=int(len(sub)), shares=sh.tolist()))
    print(f"  feeder {j + 1}: {len(sub)} labels, shares " + ", ".join(f"{s:.4f}" for s in sh))
tab = np.array([np.bincount(x[j::FEEDERS], minlength=BELTS) for j in range(FEEDERS)], dtype=float)
c, pv, _, _ = stats.chi2_contingency(tab, correction=False)
print(f"  chi-square, feeder against belt: chi2 = {c:.2f}, df = 9, p = {pv:.3f}")
res["feeder_view"] = dict(feeders=feeder, chi2=float(c), p=float(pv),
                          divisible=bool(N % FEEDERS == 0), period_per_feeder=int(N // FEEDERS))
print(f"  each feeder's own stream is periodic with period {N // FEEDERS} labels "
      f"when every feeder admits at every step")

with open(os.path.join(OUT, "sequence_dependence.json"), "w") as f:
    json.dump(res, f, indent=1)
print("\nwrote results/sequence_dependence.json")
