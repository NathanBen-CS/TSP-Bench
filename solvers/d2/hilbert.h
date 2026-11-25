#ifndef HILBERT_H
#define HILBERT_H

/* Fill tour[] / pos[] with the cities in the order a Hilbert space-filling
   curve visits them.  Geometric instances only; O(N log N). */

void hilbertTour (void);

/* Sort all cities along the Hilbert curve on a 2^order x 2^order grid:
   cities[i] = i-th city along the curve, keys[i] = its curve position
   (keys may be NULL).  Geometric instances only. */
void hilbertSort (int order, int *cities, unsigned long long *keys);

/* Position of grid cell (x, y) along the Hilbert curve on an n x n grid,
   n a power of 2.  (0,0) -> 0, and (n-1, 0) -> n*n - 1. */
unsigned long long hilbertIndex (unsigned x, unsigned y, unsigned n);

#endif
