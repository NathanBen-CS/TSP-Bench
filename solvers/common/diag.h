#ifndef DIAG_H
#define DIAG_H

#include "boxes.h"

/* How far out of its box does a city have to look?  Returns the mean share of
   each city's nbK nearest neighbours that lie outside its own box (needs
   buildNeighbors).  *anyOutside receives the same for the nearest one only. */
double outsideShare (const BoxSet *bs, double *nearestOutside);

/* Tour edges in tour[] that join two different boxes, and their length. */
long long crossingEdges (const BoxSet *bs, long long *length);

/* Print boxes, avg_box, capped_boxes, largest_box and the numbers above as
   key=value lines, prefixed with "diag_". */
void printDiagnostics (const BoxSet *bs, const BoxStats *st);

#endif
