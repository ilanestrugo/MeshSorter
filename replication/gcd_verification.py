# -*- coding: utf-8 -*-
"""
Structural verification of the residue-class / gcd claims for the MeshSorter.

The number of independent subsystems ("families") is computed by brute force
from the model's MECHANICS: one admission draw is traced from its loading
epoch through the drop opportunities it is exposed to on each primary belt
and through the recirculations it survives, and every primary slot whose
update can depend on that draw is unioned into one component.  The resulting
component count is then compared with the gcd formulae.

Nothing in build_components() uses the algebraic argument that the proof
would use; it only encodes "these slots can depend on the same draw".

Claims under test
-----------------
single-drop, loop lengths L_1..L_m   ->  g = gcd(L_2, ..., L_m)
      (L_1 excluded: feeder 1 never blocks, so it never recirculates)
dual-drop, common length L           ->  g = gcd(L, d_1 - w, ..., d_n - w)
      with d_i = q^b_i - q^f_i  and  w the feeder width
"""

from math import gcd
from functools import reduce
import itertools
import random


# ----------------------------------------------------------------- union-find
class DSU:
    def __init__(self, size):
        self.p = list(range(size))
        self.r = [0] * size

    def find(self, x):
        p = self.p
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.r[ra] < self.r[rb]:
            ra, rb = rb, ra
        self.p[rb] = ra
        if self.r[ra] == self.r[rb]:
            self.r[ra] += 1

    def count(self):
        return len({self.find(x) for x in range(len(self.p))})


def lcm(a, b):
    return a * b // gcd(a, b)


# ------------------------------------------------------------------- mechanics
def build_components(n, Ls, qf, qb=None, omega=0, dual=False,
                     feeder1_recirculates=False, qb_per_feeder=None,
                     feeder1_reaches_backward=False):
    """Count independent families by tracing admission draws.

    The draw, not the item and not the physical slot: one item carries one
    destination and so attempts one primary belt, while the draw is the state
    a feeder carries from one attempt to the next.  A slot on feeder 1 does
    come round again after L_1 steps; what does not survive the pass is the
    draw, since the item always transfers and the slot is reloaded
    independently.  That is why feeder 1 contributes no multiple of L_1.

    n      number of primary belts
    Ls     loop length of each feeder, Ls[0] is feeder 1
    qf     forward drop-point loop position on belt i (common to all feeders)
    qb     backward drop-point loop position on belt i (dual-drop only)
    omega  feeder width: a primary slot at the forward crossing at time t
           reaches the backward crossing at time t + omega

    Feeder 1 transfers at its first forward opportunity, so it neither
    recirculates nor ever reaches a backward drop point; both exclusions are
    overridable so that the controls can show they are doing real work.

    qb_per_feeder, when given, supplies a separate backward profile for each
    feeder, which is the relaxed form in which the feeders need not share a
    shape.  It is what distinguishes the corrected generating set, which runs
    over feeders 2..m, from one that also admits feeder 1's offsets.
    """
    m = len(Ls)
    M = reduce(lcm, Ls)

    # vertices: primary slot (belt i, label k mod M), k in forward-crossing time
    def v(i, k):
        return i * M + (k % M)

    dsu = DSU(n * M)

    for j in range(m):
        L = Ls[j]
        # feeder 1 never blocks, so its item makes exactly one pass and never
        # travels on from a forward drop point to the backward one
        recirculates = feeder1_recirculates or j > 0
        backward = dual and (feeder1_reaches_backward or j > 0)
        qbj = (qb_per_feeder[j] if qb_per_feeder is not None else qb)

        for e in range(M):            # loading epoch of one admission draw
            # (a) what this draw is exposed to on a single pass through the mesh
            touched = []
            for i in range(n):
                touched.append(v(i, e + qf[i]))
                if backward:
                    # the slot at the backward crossing at time e + qb[i] was
                    # at the forward crossing omega steps earlier
                    touched.append(v(i, e + qbj[i] - omega))
            first = touched[0]
            for x in touched[1:]:
                dsu.union(first, x)

            # (b) if the item can be blocked the loop slot comes round again
            #     one loop later still holding it, which links what it faces
            #     now to what it faces on the next revolution
            if recirculates:
                for i in range(n):
                    dsu.union(v(i, e + qf[i]), v(i, e + qf[i] + L))
                    if backward:
                        dsu.union(v(i, e + qbj[i] - omega),
                                  v(i, e + qbj[i] - omega + L))

    return dsu.count()


# ------------------------------------------------------------------ predictions
def predict_single(Ls):
    return reduce(gcd, Ls[1:]) if len(Ls) > 1 else 0


def predict_dual(L, qf, qb, omega):
    g = L
    for a, b in zip(qf, qb):
        g = gcd(g, abs((b - a) - omega))
    return g


# ----------------------------------------------------------------------- sweeps
def sweep_single(verbose=False):
    print('=== single-drop, unequal loop lengths ===')
    print('claim: families = gcd(L_2, ..., L_m)')
    cases = mism = 0
    rng = random.Random(20260915)

    profiles = []
    for m in (2, 3, 4):
        for _ in range(60):
            profiles.append([rng.randint(4, 14) for _ in range(m)])
    # structured profiles of the kind the paper uses
    profiles += [[6, 6, 6, 6], [6, 6, 7, 8], [6, 6, 8, 10], [6, 6, 9, 12],
                 [8, 8, 8], [8, 8, 9], [8, 8, 10], [8, 8, 12],
                 [12, 12, 13, 14], [12, 12, 14, 16], [10, 10, 12, 14]]

    for Ls in profiles:
        for n in (2, 3, 4):
            if reduce(lcm, Ls) * n > 400000:
                continue
            qf = [2 * i + 1 for i in range(n)]
            got = build_components(n, Ls, qf)
            want = predict_single(Ls)
            cases += 1
            if got != want:
                mism += 1
                print('  MISMATCH n=%d L=%s : got %d, gcd predicts %d'
                      % (n, Ls, got, want))
            elif verbose:
                print('  ok n=%d L=%-22s families=%d' % (n, Ls, got))
    print('  %d cases, %d mismatches' % (cases, mism))
    return mism


def sweep_dual(verbose=False):
    print()
    print('=== dual-drop, common loop length ===')
    print('claim: families = gcd(L, d_1 - w, ..., d_n - w)')
    cases = mism = 0
    rng = random.Random(4242)

    for _ in range(400):
        n = rng.randint(2, 4)
        m = rng.randint(2, 4)
        L = rng.randint(4 * n + 4, 4 * n + 20)
        omega = rng.randint(1, 4)
        # outbound crossings in belt order, return crossings in reverse order,
        # as in the paper's layout
        spacing = rng.randint(1, 3)
        base = rng.randint(0, 2)
        qf = [base + spacing * i for i in range(n)]
        qb = [L - 1 - base - spacing * i for i in range(n)]
        if max(qf) >= min(qb):
            continue
        Ls = [L] * m
        got = build_components(n, Ls, qf, qb, omega, dual=True)
        want = predict_dual(L, qf, qb, omega)
        cases += 1
        if got != want:
            mism += 1
            print('  MISMATCH n=%d m=%d L=%d w=%d qf=%s qb=%s : got %d, gcd predicts %d'
                  % (n, m, L, omega, qf, qb, got, want))
        elif verbose:
            print('  ok n=%d L=%-3d w=%d families=%d' % (n, L, omega, got))
    print('  %d cases, %d mismatches' % (cases, mism))
    return mism


def sweep_dual_arbitrary():
    """Same claim, but with drop points placed arbitrarily rather than in the
    paper's ordered layout, to widen coverage."""
    print()
    print('=== dual-drop, common L, arbitrary drop-point positions ===')
    cases = mism = 0
    rng = random.Random(777)
    for _ in range(400):
        n = rng.randint(2, 4)
        L = rng.randint(8, 26)
        omega = rng.randint(0, 5)
        pos = rng.sample(range(L), 2 * n)
        qf, qb = pos[:n], pos[n:]
        got = build_components(n, [L] * rng.randint(2, 3), qf, qb, omega, dual=True)
        want = predict_dual(L, qf, qb, omega)
        cases += 1
        if got != want:
            mism += 1
            print('  MISMATCH n=%d L=%d w=%d qf=%s qb=%s : got %d, predicts %d'
                  % (n, L, omega, qf, qb, got, want))
    print('  %d cases, %d mismatches' % (cases, mism))
    return mism


def sweep_dual_unequal_shapes():
    """The case that tells the two generating sets apart.

    When every feeder has the same shape, feeder 1's offsets equal those of
    feeders 2..m and including them changes nothing, which is why a sweep over
    common shapes cannot detect the error.  Give feeder 1 a shape of its own
    and the two predictions separate.  The claim is that the enumeration
    follows the one that runs over feeders 2..m only.
    """
    print()
    print('=== dual-drop, feeder 1 shaped differently from feeders 2..m ===')
    print('claim: families = gcd(L, {d_ij - w : j >= 2}), feeder 1 excluded')
    cases = mism = discriminating = 0
    rng = random.Random(4242)
    for _ in range(600):
        n = rng.randint(2, 4)
        L = rng.randint(10, 30)
        omega = rng.randint(0, 4)
        m = rng.randint(2, 4)

        def shape():
            pos = sorted(rng.sample(range(L), 2 * n))
            return pos[:n], pos[n:]

        qf, qb_rest = shape()
        _, qb_one = shape()                 # feeder 1's own backward profile
        qbs = [qb_one] + [qb_rest] * (m - 1)

        got = build_components(n, [L] * m, qf, omega=omega, dual=True,
                               qb_per_feeder=qbs)

        def predict(profiles):
            g = L
            for prof in profiles:
                for a, b in zip(qf, prof):
                    g = gcd(g, abs((b - a) - omega))
            return g

        want = predict([qb_rest])                    # feeders 2..m only
        naive = predict([qb_one, qb_rest])           # the old, wider set
        cases += 1
        if want != naive:
            discriminating += 1
        if got != want:
            mism += 1
            if mism <= 5:
                print('  MISMATCH n=%d m=%d L=%d w=%d : got %d, predicts %d'
                      % (n, m, L, omega, got, want))
    print('  %d cases, %d of them discriminating, %d mismatches'
          % (cases, discriminating, mism))
    return mism


def control_feeder1():
    """Control: if feeder 1 were allowed to recirculate, L_1 should enter the
    gcd. Confirms the exclusion of L_1 is doing real work rather than being
    vacuous."""
    print()
    print('=== control: feeder 1 forced to recirculate ===')
    checked = agree_all = agree_tail = 0
    for Ls in ([6, 9, 12], [8, 12, 20], [10, 15, 25], [4, 6, 10], [9, 6, 12]):
        n = 3
        qf = [2 * i + 1 for i in range(n)]
        got = build_components(n, Ls, qf, feeder1_recirculates=True)
        with_all = reduce(gcd, Ls)
        tail_only = reduce(gcd, Ls[1:])
        checked += 1
        agree_all += (got == with_all)
        agree_tail += (got == tail_only)
        print('  L=%-14s families=%-3d gcd(all)=%-3d gcd(L_2..L_m)=%d'
              % (Ls, got, with_all, tail_only))
    print('  matches gcd over all feeders: %d/%d ; matches tail-only: %d/%d'
          % (agree_all, checked, agree_tail, checked))




def headline_cases():
    
    print('=== consistency with Section 4.1: common length gives L families ===')
    for L in (8, 12, 24):
        for m in (3, 4):
            got = build_components(4, [L] * m, [2 * i + 1 for i in range(4)])
            print('  n=4 m=%d L=%-3d -> families=%-3d  (Section 4.1 says %d)%s'
                  % (m, L, got, L, '' if got == L else '   <-- MISMATCH'))
    
    print()
    print('=== discriminating control: does excluding L_1 change the answer? ===')
    print('    profiles chosen so gcd(L_1..L_m) != gcd(L_2..L_m)')
    for Ls in ([9, 6, 12], [4, 6, 9], [5, 10, 15], [7, 14, 21], [8, 6, 10]):
        qf = [2 * i + 1 for i in range(3)]
        real = build_components(3, Ls, qf, feeder1_recirculates=False)
        forced = build_components(3, Ls, qf, feeder1_recirculates=True)
        tail = reduce(gcd, Ls[1:])
        allg = reduce(gcd, Ls)
        ok = 'ok' if (real == tail and forced == allg) else 'MISMATCH'
        print('  L=%-13s feeder1 passive -> %-3d (gcd tail=%-3d) | forced to '
              'recirculate -> %-3d (gcd all=%-3d)  %s'
              % (Ls, real, tail, forced, allg, ok))
    
    print()
    # The enumeration builds a graph on n * lcm(L_1..L_m) vertices, so a fully
    # staggered profile at the physical length is out of reach: lcm of
    # 104,105,106,107 is 6.2e7.  The claim is structural and scale free, so the
    # staggered profiles below are scaled down, while the common-length case
    # and the dual-drop block that follows run at the physical L = 20n+24.
    print('=== single-drop profiles, 4 primary belts ===')
    for Ls in ([104, 104, 104, 104], [104, 104, 106, 108],
               [24, 24, 26, 28], [24, 24, 25, 26], [24, 25, 26, 27]):
        qf = [1 + 10 * i for i in range(4)] if Ls[0] == 104 \
            else [1 + 2 * i for i in range(4)]
        got = build_components(4, Ls, qf)
        print('  L=%-26s families=%-6d  gcd(L_2..L_4)=%d'
              % (Ls, got, predict_single(Ls)))

    print()
    print('=== dual-drop at the layout of Figure 3: L=104, n=4, w=10 ===')
    print('    qf = 1,11,21,31 and qb = 93,83,73,63, so d_i = 92,72,52,32')
    L, n = 104, 4
    for omega in (10, 8, 6, 2):
        qf = [1 + 10 * i for i in range(n)]
        for d in ([92, 72, 52, 32], [96, 76, 56, 36], [97, 77, 57, 37]):
            qb = [qf[i] + d[i] for i in range(n)]
            if max(qb) >= L:
                continue
            got = build_components(n, [L] * 4, qf, qb, omega, dual=True)
            print('  w=%-2d d=%-20s d_i-w=%-20s families=%-3d  gcd=%d'
                  % (omega, d, [x - omega for x in d], got,
                     predict_dual(L, qf, qb, omega)))


if __name__ == '__main__':
    bad = 0
    bad += sweep_single()
    bad += sweep_dual()
    bad += sweep_dual_arbitrary()
    bad += sweep_dual_unequal_shapes()
    control_feeder1()
    print()
    headline_cases()
    print()
    print('TOTAL MISMATCHES:', bad)
