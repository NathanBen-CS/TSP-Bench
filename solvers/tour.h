#ifndef TOUR_H
#define TOUR_H

#include "Instance.h"

/*********************************************************************

Shared tour representation for the solvers in this directory.

	tour[i]    city at position i of the tour, i = 0..ncities-1
	pos[c]     position of city c in tour[], c = 1..ncities

The tour is a cycle: after tour[ncities-1] comes tour[0].

 *********************************************************************/

extern int *tour;
extern int *pos;

static inline int succ (int c) {             /* city after c */
  int i = pos[c] + 1;
  return tour[i == ncities ? 0 : i];
}

static inline int pred (int c) {             /* city before c */
  int i = pos[c] - 1;
  return tour[i < 0 ? ncities - 1 : i];
}

/* How many cities apart c1 and c2 are along the tour (the shorter way round). */
static inline int tourGap (int c1, int c2) {
  int d = pos[c2] - pos[c1];
  if (d < 0) d = -d;
  return (2*d > ncities) ? ncities - d : d;
}

/* How many cities exchange2(x1,x2,y1,y2) would have to reverse (it always
   reverses the shorter side).  Local search uses this to skip moves that
   would be very slow on huge tours. */
static inline int exchangeCost (int x1, int x2, int y1, int y2) {
  int forward = (succ(x1) == x2);              /* same test as exchange2 */
  int from = forward ? x2 : x1;
  int to   = forward ? y1 : y2;
  int len  = pos[to] - pos[from];
  if (len < 0) len += ncities;
  len++;                                       /* cities on the path */
  return (2*len > ncities) ? ncities - len : len;
}

void   tourAlloc (void);
void   tourSetPositions (void);              /* rebuild pos[] from tour[] */
void   reversePath (int from, int to);
void   exchange2 (int x1, int x2, int y1, int y2);
double tourLength (void);
void   writeTour (char *filename);

#endif
