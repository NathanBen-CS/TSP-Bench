#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "instance3d.h"
#include "Utilities.h"

/*********************************************************************

Reader for 3D instances: TSPLIB files with EDGE_WEIGHT_TYPE : EUC_3D.

The 3D counterpart of codes/Instance.c.  It defines the same globals
(codes/Instance.h), so everything that only uses ncities and distance()
links with either reader, and adds zz[].

distance() is TSPLIB's nearest-integer rule, the same expression as
Concorde's euclid3d code that the Waterloo star tours were measured with
(https://math.uwaterloo.ca/tsp/star/euclid3d.txt):

	(int) (sqrt (dx*dx + dy*dy + dz*dz) + 0.5)

Differences from codes/Instance.c: coordinates are always read as real
numbers (the star files are not integral), header lines may come in any
order, unknown header lines are ignored, and the EOF line is optional
(star10k.tsp has none).

 *********************************************************************/

char *filename;
int  ncities;
char *instanceName;
int  (*distance)(int, int);
int  isEuclidean;
int  notIntegral;
double *xx, *yy, *zz;
int  **distMat;                  /* never used in 3D; declared by Instance.h */

#define LINE 512
static char name[LINE];

int distance3d (int a, int b) {
  double dx = xx[a] - xx[b];
  double dy = yy[a] - yy[b];
  double dz = zz[a] - zz[b];
  return (int) (sqrt (dx*dx + dy*dy + dz*dz) + 0.5);
}

static int isBlank (char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; }

/* If the line is "KEY : value", return the value (trimmed), else 0. */
static char *match (char *line, const char *key) {
  size_t n = strlen (key);
  char *s = line, *e;
  while (*s && isBlank (*s)) s++;
  if (strncmp (s, key, n)) return 0;
  s += n;
  if (*s && !isBlank (*s) && *s != ':') return 0;  /* a longer word */
  while (*s && isBlank (*s)) s++;
  if (*s == ':') s++;
  while (*s && isBlank (*s)) s++;
  for (e = s + strlen (s); e > s && isBlank (e[-1]); ) *--e = 0;
  return s;
}

static void readCoords (FILE *fp, long firstLine) {

  char line[LINE], msg[LINE], *p, *q;
  int i, c, k;
  double v[3];
  char *seen;
  long n = firstLine;

  xx = (double *) calloc (ncities + 1, sizeof(double));
  yy = (double *) calloc (ncities + 1, sizeof(double));
  zz = (double *) calloc (ncities + 1, sizeof(double));
  seen = (char *) calloc (ncities + 1, 1);
  if (!xx || !yy || !zz || !seen) fatal ("Out of memory reading coordinates");
  notIntegral = 0;

  for (i = 0; i < ncities; ) {
    if (!fgets (line, LINE, fp)) fatal ("File ends before all cities were read");
    n++;
    c = (int) strtol (line, &p, 10);
    if (p == line) {                               /* no number: blank line, */
      for (p = line; *p && isBlank (*p); ++p) ;    /*   or trouble           */
      if (!*p) continue;
      snprintf (msg, sizeof msg, "Line %ld: expected a city", n);
      fatal (msg);
    }
    if (c < 1 || c > ncities) fatal ("City number is out of bounds");
    if (seen[c]) fatal ("City listed twice");
    seen[c] = 1;
    for (k = 0; k < 3; ++k) {
      v[k] = strtod (p, &q);
      if (q == p) {
	snprintf (msg, sizeof msg, "Line %ld: city %d needs three coordinates", n, c);
	fatal (msg);
      }
      if (v[k] != floor (v[k])) notIntegral = 1;
      p = q;
    }
    xx[c] = v[0];
    yy[c] = v[1];
    zz[c] = v[2];
    i++;
  }
  free (seen);
  isEuclidean = 1;
  distance = distance3d;
}

void readInstance (char *s) {

  char line[LINE], type[LINE] = "", weight[LINE] = "", *v;
  long n = 0;
  FILE *fp = fopen (s, "rb");

  if (!fp) fatal ("Couldn't open input file.\n");
  filename = s;
  instanceName = name;
  strcpy (name, "NONE");
  ncities = 0;
  distance = 0;

  while (fgets (line, LINE, fp)) {
    n++;
    if      ((v = match (line, "NAME")))             strcpy (name, v);
    else if ((v = match (line, "TYPE")))             strcpy (type, v);
    else if ((v = match (line, "DIMENSION")))        ncities = atoi (v);
    else if ((v = match (line, "EDGE_WEIGHT_TYPE"))) strcpy (weight, v);
    else if (match (line, "NODE_COORD_SECTION")) {
      if (strcmp (weight, "EUC_3D")) fatal ("The 3D reader only supports EDGE_WEIGHT_TYPE EUC_3D");
      if (ncities < 1) fatal ("No DIMENSION before the coordinates");
      readCoords (fp, n);
      break;
    }
    else if (match (line, "EOF")) break;
    /* anything else (COMMENT, NODE_COORD_TYPE, ...) is ignored */
  }
  fclose (fp);
  if (strncmp (type, "TSP", 3)) fatal ("Input file is not of TYPE TSP");
  if (!distance) fatal ("No coordinates in input file");
}
