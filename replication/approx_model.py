#!/usr/bin/env python3
"""
The analytical approximation model of Section 6 of the manuscript.

The model replaces the simulation by a recursion along the feeders of a single
primary belt.  Items are assumed to arrive at each drop point independently of
one another, so that each buffer behaves as a Geo/Geo/1/c queue whose birth and
death probabilities are read off the belt utilization upstream of it.  Under
that assumption the utilization of a primary belt after the last feeder can be
computed in O(m) arithmetic operations, against the millions of time steps a
replication needs.

Everything here is a direct transcription of the equations in the manuscript.
The section numbers in the comments refer to it.

    single-drop          u_j = u_{j-1} + (1 - u_{j-1}) P_j
    dual-drop            u_jf = u_{j-1,b} + (1 - u_{j-1,b}) P_jf
                         u_jb = u_jf      + (1 - u_jf)      P_jb

with the drop probabilities P as given below.  Feeder 1 is never blocked, so
u_0 = 0 and any capacity placed at feeder 1 is ignored.

The throughput of the whole sorter is n * u_m * r; we take the belt speed r as
one slot per time step throughout, so throughput and n * u_m coincide and the
relative gap against the simulation is the same whether it is measured on
throughput or on utilization.

One correction to the manuscript is applied here.  In the case with buffers at
both drop points of the same feeder, the drop probability of the backward point
must carry the same factor pi_c of the forward buffer that its own birth
probability carries: an item reaches the backward point only if the forward
buffer could not absorb it.  Written without that factor the model violates
flow conservation, the identity sum_j w_j = n u_m that the manuscript states as
an internal consistency check, by as much as half an item per time step.  With
it the identity holds exactly in all four cases.  The published results are
unaffected, since they place buffers only at backward drop points.

Run this file directly to execute its self-checks.
"""

from __future__ import annotations

__all__ = ["utilization", "throughput", "loader_utilizations", "trace"]


# --------------------------------------------------------------------------
#  The Geo/Geo/1/c buffer
# --------------------------------------------------------------------------
def geo_geo_pi(rho: float, c: int) -> tuple[float, float]:
    """Return (pi_0, pi_c) for a Geo/Geo/1/c buffer at relative load rho.

    Equations (pi close formula) and (pi close formula full).  A buffer of
    capacity zero is always both empty and full, which is what makes the
    no-buffer case a special case of the buffered one rather than a branch.
    """
    if c <= 0:
        return 1.0, 1.0
    if rho <= 0.0:
        return 1.0, 0.0
    if abs(rho - 1.0) < 1e-12:            # the limiting uniform distribution
        return 1.0 / (c + 1), 1.0 / (c + 1)
    denom = 1.0 - rho ** (c + 1)
    if abs(denom) < 1e-300:
        return (0.0, 1.0) if rho > 1.0 else (1.0, 0.0)
    pi0 = (1.0 - rho) / denom
    pic = (1.0 - rho) * rho ** c / denom
    return min(max(pi0, 0.0), 1.0), min(max(pic, 0.0), 1.0)


# --------------------------------------------------------------------------
#  The recursion
# --------------------------------------------------------------------------
def trace(n: int, m: int, c_forward, c_backward=None):
    """Run the recursion and return the per-feeder quantities.

    n            primary belts
    m            feeder loops
    c_forward    capacity at the forward drop point of each feeder, length m
    c_backward   the same at the backward drop points; None for a single-drop
                 system, in which there is only one drop point per crossing

    Returns a dict with the utilization after each drop point and the loader
    utilization of each feeder.  Entry j-1 of each list refers to feeder j.
    """
    dual = c_backward is not None
    cf = [int(x) for x in c_forward] + [0] * (m - len(c_forward))
    cb = ([int(x) for x in c_backward] + [0] * (m - len(c_backward))) if dual \
        else [0] * m
    cf[0] = cb[0] = 0                      # feeder 1 is never blocked

    u_f, u_b, w = [], [], []
    u_prev = 0.0                           # utilization arriving at feeder j

    for j in range(m):
        # ---- forward drop point -------------------------------------------
        # The buffer there sees a birth when an item for this belt arrives
        # (probability 1/n) and the slot is taken (probability u_prev), and a
        # death when a free slot arrives and no new item does.
        rho_f = (u_prev / ((n - 1.0) * (1.0 - u_prev))) if u_prev < 1.0 else 1e9
        pi0_f, pic_f = geo_geo_pi(rho_f, cf[j])
        P_f = (1.0 - pi0_f) + pi0_f / n
        uf = u_prev + (1.0 - u_prev) * P_f

        if not dual:
            # Loader utilization, equation (workload).  With no buffer pi_c is
            # one, and the expression collapses to 1 - u_prev.
            w.append(1.0 - pic_f * u_prev)
            u_f.append(uf)
            u_prev = uf
            continue

        # ---- backward drop point ------------------------------------------
        # Only items that failed to drop forward reach it, so the arrival
        # probability carries the factor u_prev.
        # Probability that an item is waiting at the backward drop point at
        # all: it is destined for this belt (1/n), it found the slot taken at
        # the forward drop point (u_prev), and, when there is a forward buffer,
        # that buffer was full so it could not be absorbed there (pi_c).
        arriving = (u_prev / n) * (pic_f if cf[j] > 0 else 1.0)

        if cb[j] <= 0:
            P_b = arriving
            pi0_b, pic_b = 1.0, 1.0
        else:
            if cf[j] > 0:
                # Buffers at both drop points.
                b = pic_f * u_prev * uf / n
                d = (1.0 - uf) * (1.0 - u_prev * pic_f / n)
            else:
                # The configuration the design rule selects: backward only.
                b = uf * u_prev / n
                d = (1.0 - uf) * ((1.0 - 1.0 / n) + (1.0 - u_prev) / n)
            rho_b = (b / d) if d > 0.0 else 1e9
            pi0_b, pic_b = geo_geo_pi(rho_b, cb[j])
            P_b = (1.0 - pi0_b) + pi0_b * arriving

        ub = uf + (1.0 - uf) * P_b
        # Equation (workload dual).  An item is turned away only if the slot
        # is taken at the forward drop point AND the forward buffer, if any,
        # is full AND the slot is still taken at the backward drop point AND
        # the backward buffer is full.  With no forward buffer pi_c is one and
        # the product collapses to the form printed in the manuscript.
        w.append(1.0 - u_prev * pic_f * uf * pic_b)
        u_f.append(uf)
        u_b.append(ub)
        u_prev = ub

    return {"u_forward": u_f, "u_backward": u_b if dual else None,
            "u_final": (u_b[-1] if dual else u_f[-1]), "loaders": w}


def utilization(n: int, m: int, c_forward, c_backward=None) -> float:
    """Utilization of a primary belt as it leaves the feeding area."""
    return trace(n, m, c_forward, c_backward)["u_final"]


def throughput(n: int, m: int, c_forward, c_backward=None) -> float:
    """Items per time step delivered by the whole sorter, at belt speed one."""
    return n * utilization(n, m, c_forward, c_backward)


def loader_utilizations(n: int, m: int, c_forward, c_backward=None):
    """Utilization of each loading station, the analytical counterpart of the
    measurements in Section 5.3."""
    return trace(n, m, c_forward, c_backward)["loaders"]


# --------------------------------------------------------------------------
#  Self-checks
# --------------------------------------------------------------------------
def _selfcheck() -> None:
    import itertools

    # 1.  Single-drop without buffers has the closed form 1 - (1 - 1/n)^m.
    for n in range(2, 10):
        for m in range(1, 10):
            got = utilization(n, m, [0] * m)
            want = 1.0 - (1.0 - 1.0 / n) ** m
            assert abs(got - want) < 1e-12, (n, m, got, want)

    # 2.  The manuscript states that the loader utilizations sum to n * u_m
    #     exactly, in every configuration.  This is the model's own internal
    #     consistency check, and it exercises every one of the four dual-drop
    #     cases.
    worst, seen = 0.0, {"single": 0, "none": 0, "forward": 0,
                        "backward": 0, "both": 0}
    for n in (3, 4, 6):
        for m in (3, 4):
            for cf in itertools.product(range(3), repeat=m - 1):
                cfv = [0] + list(cf)
                t = trace(n, m, cfv)
                worst = max(worst, abs(sum(t["loaders"]) - n * t["u_final"]))
                seen["single"] += 1
                for cb in itertools.product(range(3), repeat=m - 1):
                    cbv = [0] + list(cb)
                    t = trace(n, m, cfv, cbv)
                    worst = max(worst, abs(sum(t["loaders"]) - n * t["u_final"]))
                    f, b = max(cfv) > 0, max(cbv) > 0
                    seen["both" if f and b else "forward" if f
                         else "backward" if b else "none"] += 1
    assert worst < 1e-9, worst
    assert all(v > 0 for v in seen.values()), seen

    # 3.  Adding capacity never lowers the predicted utilization, and moving
    #     capacity from a forward drop point to the backward one never lowers
    #     it either.  These are the two empirical design rules of Section 5.1,
    #     which the model should reproduce rather than contradict.
    for n in (4,):
        for m in (3, 4, 5):
            for c in range(4):
                base = [0] + [c] * (m - 1)
                more = [0] + [c + 1] * (m - 1)
                assert utilization(n, m, base) <= utilization(n, m, more) + 1e-12
                fwd = utilization(n, m, [0] + [c] * (m - 1), [0] * m)
                bwd = utilization(n, m, [0] * m, [0] + [c] * (m - 1))
                assert bwd >= fwd - 1e-12, (n, m, c, fwd, bwd)

    print("approx_model: self-checks passed on "
          f"{sum(seen.values())} configurations "
          f"({', '.join(f'{k} {v}' for k, v in seen.items())}); "
          f"largest flow-conservation residual {worst:.2e}")


if __name__ == "__main__":
    _selfcheck()
