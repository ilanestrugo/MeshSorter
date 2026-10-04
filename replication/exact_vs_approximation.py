# -*- coding: utf-8 -*-
"""Is the approximation an upper bound on the exact throughput?

Exact Q(n,m): stationary law of the lumped partition chain of Section 4.3.
Approximation A(n,m) = n(1-(1-1/n)^m), which is also the expected number of
occupied bins when m balls are thrown independently into n bins.

Also checks the per-feeder reduction:  Q <= A  for all m  <=>  w_j <= (1-1/n)^(j-1)
for all j, where w_j = Q(n,j) - Q(n,j-1) is the loader utilization of feeder j.
"""
import sys
from functools import lru_cache


def partitions(m, maxparts):
    def rec(rem, cap, parts):
        if rem == 0:
            yield tuple(parts); return
        if len(parts) == maxparts:
            return
        for v in range(min(rem, cap), 0, -1):
            yield from rec(rem - v, v, parts + [v])
    return list(rec(m, m, []))


def transition(lam, n):
    """One row, built by taking the K transferring feeders one at a time.
    State during the build: (block sizes on used destinations, #used)."""
    K = len(lam)
    surv = tuple(sorted([b - 1 for b in lam if b >= 2], reverse=True))
    cur = {(surv, len(surv)): 1.0}
    for _ in range(K):
        nxt = {}
        for (blocks, used), p in cur.items():
            for i in range(len(blocks)):          # join an existing block
                nb = list(blocks); nb[i] += 1
                k = (tuple(sorted(nb, reverse=True)), used)
                nxt[k] = nxt.get(k, 0.0) + p / n
            if used < n:                          # open a free destination
                nb = tuple(sorted(list(blocks) + [1], reverse=True))
                k = (nb, used + 1)
                nxt[k] = nxt.get(k, 0.0) + p * (n - used) / n
        cur = nxt
    row = {}
    for (blocks, _used), p in cur.items():
        k = tuple(b for b in blocks if b > 0)
        row[k] = row.get(k, 0.0) + p
    return row


@lru_cache(maxsize=None)
def Q(n, m):
    if m == 0:
        return 0.0
    states = partitions(m, n)
    idx = {s: i for i, s in enumerate(states)}
    N = len(states)
    # power iteration on the lumped chain (it has a unique closed class)
    P = [transition(s, n) for s in states]
    pi = [1.0 / N] * N
    for _ in range(20000):
        new = [0.0] * N
        for i, row in enumerate(P):
            pv = pi[i]
            if pv:
                for k, p in row.items():
                    new[idx[k]] += pv * p
        s = sum(new)
        new = [x / s for x in new]
        if max(abs(a - b) for a, b in zip(new, pi)) < 1e-15:
            pi = new; break
        pi = new
    return sum(pi[idx[s]] * len(s) for s in states)


def A(n, m):
    return n * (1.0 - (1.0 - 1.0 / n) ** m)


if __name__ == '__main__':
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12

    print('=== validate against Table 2(a) ===')
    ok = True
    for (n, m, want) in [(4, 4, 2.6210), (4, 9, 3.3444), (9, 4, 3.3440),
                         (6, 6, 3.7809), (9, 9, 5.5312), (5, 7, 3.6486)]:
        got = Q(n, m)
        good = abs(got - want) < 6e-4
        ok &= good
        print('  Q(%d,%d) = %.4f   table %.4f   %s'
              % (n, m, got, want, 'ok' if good else 'MISMATCH'))
    print('  validation:', 'PASSED' if ok else 'FAILED')

    print()
    print('=== per-feeder reduction:  w_j  vs  (1-1/n)^(j-1) ===')
    for n in (4, 6, 9):
        print('  n=%d' % n)
        for j in range(1, 8):
            w = Q(n, j) - Q(n, j - 1)
            a = (1.0 - 1.0 / n) ** (j - 1)
            if abs(w - a) < 1e-12:
                flag = 'EQUAL'
            elif w < a:
                flag = 'ok  (gap %.4f)' % (a - w)
            else:
                flag = '*** VIOLATION ***'
            print('    j=%d  w_j=%.6f  (1-1/n)^(j-1)=%.6f   %s' % (j, w, a, flag))

    print()
    print('=== Q(n,m) <= A(n,m) over a grid ===')
    viol, cells, worst = 0, 0, (0.0, None)
    for n in range(2, NMAX + 1):
        for m in range(1, NMAX + 1):
            q, a = Q(n, m), A(n, m)
            cells += 1
            if q > a + 1e-12:
                viol += 1
                print('  VIOLATION n=%d m=%d: Q=%.6f > A=%.6f' % (n, m, q, a))
            rel = (a - q) / q * 100
            if rel > worst[0]:
                worst = (rel, (n, m))
    print('  %d cells, %d violations; largest overshoot %.2f%% at n=%d m=%d'
          % (cells, viol, worst[0], worst[1][0], worst[1][1]))
