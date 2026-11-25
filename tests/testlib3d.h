#ifndef TESTLIB3D_H
#define TESTLIB3D_H

/* Shared by unit3d.c and boxtest3d.c: reporting, and 3D instances made in
   memory with the bench's rounding (EUC_3D, nearest integer). */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "instance3d.h"
#include "Utilities.h"

static int failures = 0;

static void report (const char *name, int ok, const char *detail) {
  printf ("%-4s  %s%s%s\n", ok ? "PASS" : "FAIL", name, detail[0] ? " - " : "", detail);
  if (!ok) failures++;
}

static unsigned long long rs = 88172645463325252ULL;          /* xorshift64 */
static unsigned long long rnd (void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static double urand (void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }
static double nrand (void) {                                   /* standard normal */
  return sqrt (-2 * log (urand() + 1e-300)) * cos (6.283185307179586 * urand());
}

enum { UNIFORM, CLUSTERED, DUPLICATES, NTYPES };
static const char *typeName[] = { "uniform", "clustered", "duplicates" };

/* Make an n-city instance in the globals instance3d.c normally fills. */
static void makeInstance (int n, int type) {

  int c, k, nc;
  double *cx, *cy, *cz, sd;

  free (xx); free (yy); free (zz);
  ncities = n;
  xx = (double *) calloc (n + 1, sizeof(double));
  yy = (double *) calloc (n + 1, sizeof(double));
  zz = (double *) calloc (n + 1, sizeof(double));
  isEuclidean = 1;
  distance = distance3d;

  if (type == UNIFORM) {                       /* like portgen3d */
    for (c = 1; c <= n; ++c) {
      xx[c] = rnd() % 1000000; yy[c] = rnd() % 1000000; zz[c] = rnd() % 1000000;
    }
  }
  else if (type == CLUSTERED) {                /* like portcgen3d: n/100 centres, */
    nc = n / 100 > 0 ? n / 100 : 1;            /* normal spread 1e6 / cbrt(n)     */
    cx = (double *) malloc (nc * sizeof(double));
    cy = (double *) malloc (nc * sizeof(double));
    cz = (double *) malloc (nc * sizeof(double));
    for (k = 0; k < nc; ++k) { cx[k] = urand() * 1e6; cy[k] = urand() * 1e6; cz[k] = urand() * 1e6; }
    sd = 1e6 / cbrt ((double) n);
    for (c = 1; c <= n; ++c) {
      k = rnd() % nc;
      xx[c] = (int) (cx[k] + sd * nrand());
      yy[c] = (int) (cy[k] + sd * nrand());
      zz[c] = (int) (cz[k] + sd * nrand());
    }
    free (cx); free (cy); free (cz);
  }
  else {                                       /* many exact duplicates: */
    nc = n / 25 > 1 ? n / 25 : 2;              /*   ~25 cities per spot  */
    for (c = 1; c <= n; ++c) {
      k = rnd() % nc;
      xx[c] = (k * 7919) % 1000;
      yy[c] = (k * 104729) % 1000;
      zz[c] = (k * 1299709) % 1000;
    }
  }
}

static int isPermutation (const int *a, int n) {   /* cities 1..n, each once */
  int i, ok = 1;
  char *seen = (char *) calloc (n + 1, 1);
  for (i = 0; i < n; ++i) {
    if (a[i] < 1 || a[i] > n || seen[a[i]]) { ok = 0; break; }
    seen[a[i]] = 1;
  }
  free (seen);
  return ok;
}

static long long cycleLength (const int *a, int n) {
  long long L = 0;
  int i;
  for (i = 0; i < n; ++i) L += distance (a[i], a[(i + 1) % n]);
  return L;
}

static int finish (void) {
  printf ("\n%s: %d failure(s)\n", failures ? "FAILED" : "ALL PASSED", failures);
  return failures ? 1 : 0;
}

#endif
