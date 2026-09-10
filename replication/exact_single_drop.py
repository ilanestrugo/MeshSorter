"""Exact steady-state throughput of the unbuffered single-drop MeshSorter.

The unbuffered single-drop system with a common feeder-loop length decomposes
into L independent copies of a "round chain" on the destinations currently held
by the m feeders.  In one round, feeder j hands its item to the primary belt
unless some feeder ahead of it in the round holds the same destination, in which
case it keeps the item for another round.  The throughput per time step is the
expected number of distinct destinations among the m held items in steady state.

The chain on {1..n}^m is lumpable onto the integer partitions of m, of which
there are p(m): 30 for m = 9, against n^m = 387,420,489 raw states.  This module
builds the lumped transition kernel exactly in rational arithmetic, solves for
the stationary distribution by exact Gaussian elimination, and returns a
Fraction, so the result carries no numerical error at all.

See Section 4 of the manuscript.

Usage
    python3 exact_single_drop.py            print the 6x6 grid used in Table 1(a)
    python3 exact_single_drop.py N M        print the exact value for one cell
    python3 exact_single_drop.py --verify   check small cases against a brute
                                            force chain on all n^m states
"""
from fractions import Fraction
import itertools, sys

def parts(m):
    """integer partitions of m as sorted tuples (descending)"""
    def rec(rem, mx):
        if rem == 0: yield ()
        for v in range(min(rem, mx), 0, -1):
            for t in rec(rem-v, v): yield (v,)+t
    return sorted(set(rec(m, m)))

def kernel(P, n):
    """P = block sizes.  Returns {next_partition: probability} (exact Fractions)."""
    k = len(P)                                   # blocks = successes = feeders redrawn
    surv = [b-1 for b in P if b >= 2]            # groups that keep a member (and a colour)
    s = len(surv)
    if s > n: raise ValueError
    # DP over the k redrawn feeders.
    # state: (tuple of extra counts on each survivor box, sorted tuple of free-box counts)
    start = ((0,)*s, ())
    cur = {start: Fraction(1)}
    inv = Fraction(1, n)
    for _ in range(k):
        nxt = {}
        for (sc, fc), p in cur.items():
            for i in range(s):                                   # join survivor box i
                t = list(sc); t[i] += 1
                key = (tuple(t), fc); nxt[key] = nxt.get(key, 0) + p*inv
            for i in range(len(fc)):                             # join an open free box
                t = list(fc); t[i] += 1
                key = (sc, tuple(sorted(t, reverse=True))); nxt[key] = nxt.get(key, 0) + p*inv
            free_left = n - s - len(fc)                           # open a brand-new box
            if free_left > 0:
                key = (sc, tuple(sorted(fc+(1,), reverse=True)))
                nxt[key] = nxt.get(key, 0) + p*inv*free_left
        cur = nxt
    out = {}
    for (sc, fc), p in cur.items():
        blocks = tuple(sorted([surv[i]+sc[i] for i in range(s)] + list(fc), reverse=True))
        out[blocks] = out.get(blocks, 0) + p
    return out

def exact(n, m):
    S = [P for P in parts(m) if len(P) <= n]
    idx = {P: i for i, P in enumerate(S)}
    K = [kernel(P, n) for P in S]
    N = len(S)
    # solve pi (I-K) = 0, sum pi = 1, exactly by Gaussian elimination over Q
    A = [[Fraction(0)]*N for _ in range(N)]        # A[c][r] : column c equation
    for i, row in enumerate(K):
        for P2, p in row.items(): A[idx[P2]][i] += p
    for c in range(N): A[c][c] -= 1
    A[N-1] = [Fraction(1)]*N                        # replace last equation by normalization
    b = [Fraction(0)]*N; b[N-1] = Fraction(1)
    for col in range(N):
        piv = next(r for r in range(col, N) if A[r][col] != 0)
        A[col], A[piv] = A[piv], A[col]; b[col], b[piv] = b[piv], b[col]
        inv = Fraction(1)/A[col][col]
        A[col] = [x*inv for x in A[col]]; b[col] *= inv
        for r in range(N):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [x - f*y for x, y in zip(A[r], A[col])]
                b[r] -= f*b[col]
    pi = b
    return sum(p*len(P) for p, P in zip(pi, S))

# --------------------------------------------------------------------------
#  Brute force reference: the chain on all n^m states, for verification only.
# --------------------------------------------------------------------------

def brute(n, m):
    """Same quantity, computed on the raw n^m chain.  Only usable for tiny n, m."""
    states = list(itertools.product(range(n), repeat=m))
    index = {s: k for k, s in enumerate(states)}
    N = len(states)
    inv = Fraction(1, n)
    rows = [dict() for _ in range(N)]
    for s in states:
        seen = set()
        keep = []                       # feeders that do NOT transfer this round
        for j, d in enumerate(s):
            if d in seen: keep.append(j)
            else: seen.add(d)
        redraw = [j for j in range(m) if j not in keep]
        base = list(s)
        for combo in itertools.product(range(n), repeat=len(redraw)):
            t = list(base)
            for j, v in zip(redraw, combo): t[j] = v
            k = index[tuple(t)]
            rows[index[s]][k] = rows[index[s]].get(k, 0) + inv ** len(redraw)
    # power iteration in exact arithmetic is slow; solve the linear system
    A = [[Fraction(0)] * N for _ in range(N)]
    for s_i, row in enumerate(rows):
        for k, p in row.items(): A[k][s_i] += p
    for c in range(N): A[c][c] -= 1
    A[N - 1] = [Fraction(1)] * N
    b = [Fraction(0)] * N; b[N - 1] = Fraction(1)
    for col in range(N):
        piv = next(r for r in range(col, N) if A[r][col] != 0)
        A[col], A[piv] = A[piv], A[col]; b[col], b[piv] = b[piv], b[col]
        f0 = Fraction(1) / A[col][col]
        A[col] = [x * f0 for x in A[col]]; b[col] *= f0
        for r in range(N):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[col])]
                b[r] -= f * b[col]
    return sum(p * len(set(s)) for p, s in zip(b, states))


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "--verify":
        for (n, m) in [(2, 2), (3, 3), (2, 4), (4, 3)]:
            a, b = exact(n, m), brute(n, m)
            print(f"n={n} m={m}  lumped={a}  brute={b}  {'ok' if a == b else 'MISMATCH'}")
    elif len(args) == 2:
        n, m = int(args[0]), int(args[1])
        v = exact(n, m)
        print(f"n={n} m={m}  {v}  =  {float(v):.9f}")
    else:
        print("Exact throughput of the unbuffered single-drop MeshSorter")
        print("rows: primary belts n,  columns: feeder loops m\n")
        print("     " + "".join(f"{m:>10}" for m in range(4, 10)))
        for n in range(4, 10):
            print(f"n={n} " + "".join(f"{float(exact(n, m)):>10.4f}" for m in range(4, 10)))
