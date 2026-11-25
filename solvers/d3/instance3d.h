#ifndef INSTANCE3D_H
#define INSTANCE3D_H

/* The shared reader interface (codes/Instance.h: ncities, distance(), xx, yy,
   readInstance(), ...) plus the third coordinate. */

#include "Instance.h"

extern double *zz;

/* TSPLIB EUC_3D distance: straight line, rounded to the nearest integer. */
int distance3d (int a, int b);

#endif
