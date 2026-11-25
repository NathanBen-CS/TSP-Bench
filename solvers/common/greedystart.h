#ifndef GREEDYSTART_H
#define GREEDYSTART_H

/* Fill tour[] / pos[] with a greedy (multi-fragment) tour built from the
   neighbour lists.  Needs buildNeighbors() to have been called first.
   Works for both geometric and matrix instances. */

void greedyTour (int candidatesPerCity);

#endif
