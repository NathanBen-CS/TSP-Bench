#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H

/* Improve the tour in tour[]/pos[] with 2-opt moves, plus Or-opt moves if
   useOrOpt is non-zero, until no improving move is left.  Needs
   buildNeighbors() to have been called.  Counts of applied moves are
   returned through the two pointers. */

void localSearch (int useOrOpt, long long *twoOptMoves, long long *orOptMoves);

/* Same moves, but repeated (every city re-queued) until a whole round makes
   no move, so no allowed improving move is left.  rounds = rounds run. */

void localSearchToOptimum (int useOrOpt, long long *twoOptMoves,
                           long long *orOptMoves, int *rounds);

#endif
