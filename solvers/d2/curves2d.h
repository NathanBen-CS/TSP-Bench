#ifndef CURVES2D_H
#define CURVES2D_H

#include "curves.h"

/*********************************************************************

The 2D curves (d2/curves2d.c).

	base 4  (Hilbert, Moore, Morton, Strip): a cell splits into 4 squares
	base 9  (Peano):                         a cell splits into 3 x 3 squares
	base 2  (Sierpinski, Knopp):             a triangle splits into 2 triangles

 *********************************************************************/

enum {
  CURVE_HILBERT_2D, CURVE_MOORE, CURVE_SIERPINSKI, CURVE_KNOPP,
  CURVE_PEANO, CURVE_MORTON, CURVE_STRIP, NCURVES
};

/* Key of the point (u, v) in the unit square [0,1] x [0,1].            */
unsigned long long curveKey (int curve, double u, double v);

/* ---- small-grid versions, used by the tests ---- */

/* Key of grid cell (x, y) on a side x side grid (side = base^(levels/..)):
   square-cell curves only.  levels = number of base-`base` digits.     */
unsigned long long curveGridKey (int curve, unsigned x, unsigned y, int levels);

/* Side of the square grid a square-cell curve uses for `levels` digits. */
unsigned long long curveGridSide (int curve, int levels);

/* Triangle curves: the corners of the triangle cell with key `key`
   among the 2^levels cells of depth `levels`.                          */
void curveTriangle (int curve, unsigned long long key, int levels, double tri[3][2]);

#endif
