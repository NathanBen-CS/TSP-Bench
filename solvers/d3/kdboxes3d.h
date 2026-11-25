#ifndef KDBOXES3D_H
#define KDBOXES3D_H

#include "boxes.h"

/*********************************************************************

kd-median boxes (--split kd): boxes that really hold t cities.

The curve-cell boxes of boxBuild have at most t cities, but on average far
fewer, because a cell splits into 8.  Here the cities are cut at the median
of the widest axis until the pieces hold exactly t cities each (one piece
holds the remainder).  The pieces are then put in the order in which
`curve` visits their centres, so M1-M3 have a chain to follow.

Held-Karp, merges and cleanup are unchanged.

 *********************************************************************/

void boxBuildKd (int curve, int t, BoxSet *bs, BoxStats *st);

#endif
