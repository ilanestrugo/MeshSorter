# -*- coding: utf-8 -*-
"""Three-way comparison on the one family where all three are available:
the unbuffered single-drop MeshSorter with a common loop length.

    exact          partition chain, Proposition 6 of the manuscript
    simulation     meshsorter, the protocol of Table 2
    approximation  n(1-(1-1/n)^m)

Table 4 of the manuscript.  The simulation uses the geometry of the paper, the
one in geometry.py, passed as --spacing 10 --turn 12 --df 14 --width 10, so that
the loops hold 20n + 24 slots; T = 1,330,000, R = 30, the warm-up rule of
Section S2 of the Supplemental Online Material, and seed 20260901.  By
Proposition 6 the exact value does not depend on the geometry, but the reported
sample comes from the protocol the paper states.  Results are written to
results/three_way/threeway_results.json.
"""
import json, subprocess, sys, os
from exact_vs_approximation import Q, A

BIN = os.environ.get("MESHSORTER_BIN",
                     "./meshsorter")
STEPS, REPS, SEED = 1_330_000, 30, 20260901
BELTS = FEEDERS = [3, 4, 5, 6, 7, 8, 9]
# geometry of Section 2: spacing 10, turn 12, so L = 20n + 24.
# By Proposition 6 the expectation does not depend on it, but the
# reported sample should come from the protocol the paper states.


def simulate(n, m):
    out = subprocess.run(
        [BIN, '-n', str(n), '-m', str(m), '--single',
         '--spacing', '10', '--turn', '12', '--df', '14', '--width', '10',
         '-T', str(STEPS), '-R', str(REPS), '--seed', str(SEED), '--json'],
        capture_output=True, text=True, check=True).stdout
    return json.loads(out)


rows = []
for n in BELTS:
    for m in FEEDERS:
        r = simulate(n, m)
        q, a = Q(n, m), A(n, m)
        s, hw = r['throughput'], r['halfwidth']
        rows.append(dict(n=n, m=m, exact=q, sim=s, hw=hw, approx=a,
                         sim_gap=s - q, approx_gap=a - q,
                         approx_rel=(a - q) / q * 100))
        print('n=%d m=%d  exact %.4f  sim %.4f +-%.5f  approx %.4f  '
              'sim-exact %+.5f  approx-exact %+.2f%%'
              % (n, m, q, s, hw, a, s - q, (a - q) / q * 100),
              flush=True)

_out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'three_way')
os.makedirs(_out, exist_ok=True)
with open(os.path.join(_out, 'threeway_results.json'), 'w') as f:
    json.dump(rows, f, indent=1)

print()
print('--- summary ---')
print('simulation vs exact : largest |gap| = %.5f  (mean %+.6f)'
      % (max(abs(r['sim_gap']) for r in rows),
         sum(r['sim_gap'] for r in rows) / len(rows)))
print('                      largest half-width = %.5f'
      % max(r['hw'] for r in rows))
print('                      cells where exact lies in the CI: %d of %d'
      % (sum(1 for r in rows if abs(r['sim_gap']) <= r['hw']), len(rows)))
print('approximation vs exact: %.2f%% to %.2f%%, above in %d of %d cells'
      % (min(r['approx_rel'] for r in rows), max(r['approx_rel'] for r in rows),
         sum(1 for r in rows if r['approx_gap'] > 0), len(rows)))
