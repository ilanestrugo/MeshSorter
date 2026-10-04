# -*- coding: utf-8 -*-
"""The feeding-area geometry, in slots, shared by every table script.

Everything imports these constants, so the geometry cannot drift between one
table and another.  They are the layout costed in Supplement S9 and drawn in
Figure 3 of the manuscript, converted at the pitch of 0.5 m used throughout:

    a primary belt with its bins on both sides           3 m
    a forklift aisle                                     2 m
      -> consecutive primary belts sit 5 m apart      = 10 slots  (SPACING)
    n belts and the closing aisle span 5n + 2 m
    each end of a feeder loop, the loop being 5 m across 4 m
    the closing aisle, charged to the ends               2 m
      -> each end of the loop consumes 6 m            = 12 slots  (TURN)

    loop length  L = n*DP + 2*TURN = 20n + 24 slots

    n = 4  ->  104 slots = 52 m, the loop of Figure 3
    n = 9  ->  204 slots = 102 m

A feeder loop is 5 m across, so the forward and backward crossings of one
feeder on a primary belt are 5 m = 10 slots apart: WIDTH = 10.

Consecutive feeders sit 5 m apart with a 2 m gap: DF = 14 slots.  By
Proposition 2 the steady state does not depend on DF at all; it is set
correctly for the record rather than because it changes an answer.

Earlier versions of these scripts used DP = 4, TURN = 4, WIDTH = 2, DF = 4,
which is a loop of 4n + 8 slots and puts consecutive primary belts 1 m apart.
The certification scripts used SPACING = 4, a loop of 8n + 8.  Neither is
buildable at the bin and aisle dimensions the paper costs.
"""

SPACING = 10              # slots between consecutive primary belts on one run
DP = 2 * SPACING          # loop slots consumed by each primary belt
TURN = 12                 # loop slots consumed by each end of the loop
WIDTH = 10                # forward-to-backward distance along a primary belt
DF = 14                   # spacing of consecutive feeders along a primary belt


def loop_length(n):
    """Base feeder-loop length for n primary belts, in slots."""
    return n * DP + 2 * TURN


if __name__ == '__main__':
    print('SPACING %d  DP %d  TURN %d  WIDTH %d  DF %d'
          % (SPACING, DP, TURN, WIDTH, DF))
    for n in range(3, 16):
        print('  n=%-3d L = %4d slots = %5.1f m' % (n, loop_length(n),
                                                    loop_length(n) * 0.5))
