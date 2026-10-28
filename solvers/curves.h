#ifndef CURVES_H
#define CURVES_H

/*********************************************************************

Space-filling curves for the box solver (2D).

Every curve gives each point of the plane a KEY: its position along the
curve.  Keys are written with `digits` digits in base `base`, and the
curve visits cells so that all points of one cell share a key prefix:

	base 4  (Hilbert, Moore, Morton, Strip): a cell splits into 4 squares
	base 9  (Peano):                         a cell splits into 3 x 3 squares
	base 2  (Sierpinski, Knopp):             a triangle splits into 2 triangles

So a box (a cell of the curve's own tree) is always one unbroken run of
the key-sorted city list, whatever the curve.  See curves.c.

 *********************************************************************/

typedef enum {
  CURVE_HILBERT, CURVE_MOORE, CURVE_SIERPINSKI, CURVE_KNOPP,
  CURVE_PEANO, CURVE_MORTON, CURVE_STRIP, NCURVES
} CurveId;

typedef struct {
  const char *name;
  int base;          /* children per cell: 4, 9 or 2              */
  int digits;        /* key length in base-`base` digits (depth)   */
  int closed;        /* 1 if the curve ends next to where it starts */
  int triangles;     /* 1 if its cells are triangles, 0 if squares  */
} CurveInfo;

extern const CurveInfo curveInfo[NCURVES];

int curveByName (const char *name);          /* -1 if unknown */

/* Key of the point (u, v) in the unit square [0,1] x [0,1].            */
unsigned long long curveKey (int curve, double u, double v);

/* Sort all cities by key: cities[i] = i-th city along the curve, keys[i]
   its key (keys may be NULL).  The bounding box is scaled into the unit
   square with one scale for both axes.  Geometric instances only.      */
void curveSort (int curve, int *cities, unsigned long long *keys);

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
