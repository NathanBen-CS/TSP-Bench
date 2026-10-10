#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H

/* Improve the tour in tour[]/pos[] with 2-opt moves, plus Or-opt moves if
   useOrOpt is non-zero, until no improving move is left.  Needs
   buildNeighbors() to have been called.  Counts of applied moves are
   returned through the two pointers. */

void localSearch (int useOrOpt, long long *twoOptMoves, long long *orOptMoves);

#endif
