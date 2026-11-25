#ifndef CURVES_H
#define CURVES_H

/*********************************************************************

Space-filling curves: the part every dimension shares.

A curve gives each point a KEY, its position along the curve.  Keys are
written with `digits` digits in base `base`, and the curve visits cells
so that all points of one cell share a key prefix.  So a box (a cell of
the curve's own tree) is always one unbroken run of the key-sorted city
list, whatever the curve and whatever the dimension.

The curves themselves differ per dimension.  This header has two
implementations, and the Makefile links one of them:

	d2/curves2d.c    7 curves in the plane   (extra test hooks: curves2d.h)
	d3/curves3d.c    curves in space         (extra test hooks: curves3d.h)

 *********************************************************************/

typedef struct {
  const char *name;
  int base;          /* children per cell                          */
  int digits;        /* key length in base-`base` digits (depth)   */
  int closed;        /* 1 if the curve ends next to where it starts */
  int triangles;     /* 1 if its cells are triangles / tetrahedra,
			0 if squares / cubes                        */
} CurveInfo;

extern const CurveInfo curveInfo[];
extern const int nCurves;                    /* entries in curveInfo[]  */

#define CURVE_HILBERT 0                      /* the default curve: entry 0
						in both dimensions       */

int curveByName (const char *name);          /* -1 if unknown */

/* Sort all cities by key: cities[i] = i-th city along the curve, keys[i]
   its key (keys may be NULL).  The bounding box is scaled into the unit
   square / cube with one scale for all axes.  Geometric instances only. */
void curveSort (int curve, int *cities, unsigned long long *keys);

#endif
