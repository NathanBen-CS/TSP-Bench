#ifndef CURVES3D_H
#define CURVES3D_H

#include "curves.h"

/*********************************************************************

The 3D curves (d3/curves3d.c).

	base 8  (Hilbert3, Morton3, Snake3, Moore3): a cell splits into 8 cubes
	base 27 (Peano3):                            a cell splits into 3 x 3 x 3
	base 2  (Sierpinski3):                       a tetrahedron splits in half

 *********************************************************************/

enum {
  CURVE_HILBERT3, CURVE_MORTON3, CURVE_PEANO3, CURVE_SNAKE3, CURVE_MOORE3,
  CURVE_SIERPINSKI3, NCURVES3
};

/* Key of the point (x, y, z), in the frame curveSort uses: the bounding box
   of all cities, one scale for the three axes.  Call curveFrame() first.   */
void curveFrame (void);
unsigned long long curveKeyAt (int curve, double x, double y, double z);

/* ---- small-grid versions, used by the tests ---- */

/* Key of grid cell (x, y, z) on a side^3 grid (cube-cell curves only).
   levels = number of base-`base` digits.                                */
unsigned long long curveGridKey3 (int curve, unsigned x, unsigned y, unsigned z, int levels);

/* Side of the grid a cube-cell curve uses for `levels` digits. */
unsigned long long curveGridSide3 (int curve, int levels);

/* Sierpinski3: the corners of the tetrahedron with key `key` among the cells
   at depth `levels` (levels counts the 3 digits that pick the first
   tetrahedron).  Returns 0 if there is no such cell.                      */
int curveTetra (unsigned long long key, int levels, double tet[4][3]);

/* Sierpinski3 key of a point of the unit cube with `levels` digits. */
unsigned long long curveTetraKey (const double p[3], int levels);

#endif
