#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Instance.h"
#include "Utilities.h"
#include "tour.h"
#include "neighbors.h"
#include "boxes.h"
#include "curves.h"
#include "diag.h"

/*********************************************************************

Box diagnostics for any instance, optionally against a finished tour:

	boxdiag instance_file tour_file|- [--curve NAME] [--t T] [--k K]

Splits the instance into boxes exactly as boxtsp does and prints boxes,
average cities per box, depth-capped boxes, the share of each city's K
nearest neighbours outside its own box and, with a tour file, how many
tour edges cross box borders.  The same program in 2D (boxdiag2d) and 3D
(boxdiag3d), so the two can be compared.

 *********************************************************************/

int main (int argc, char **argv) {

  int i, c, curve = CURVE_HILBERT, t = 10, k = 16;
  BoxSet bs;
  BoxStats st;
  FILE *fp;

  if (argc < 3) fatal ("Usage: boxdiag instance_file tour_file|- [--curve NAME] [--t T] [--k K]");
  for (i = 3; i + 1 < argc; i += 2) {
    if      (!strcmp (argv[i], "--t")) t = atoi (argv[i+1]);
    else if (!strcmp (argv[i], "--k")) k = atoi (argv[i+1]);
    else if (!strcmp (argv[i], "--curve")) {
      if ((curve = curveByName (argv[i+1])) < 0) fatal ("unknown curve");
    }
    else fatal ("unknown option");
  }
  if (t < 2 || t > BOX_MAXT) fatal ("--t must be between 2 and 12");

  readInstance (argv[1]);
  buildNeighbors (k);
  boxBuild (curve, t, &bs, &st);
  if (strcmp (argv[2], "-")) {
    tourAlloc();
    fp = fopen (argv[2], "r");
    if (!fp) fatal ("Can't open tour file");
    for (i = 0; i < ncities; ++i) {
      if (fscanf (fp, "%d", &c) != 1 || c < 1 || c > ncities) fatal ("bad tour file");
      tour[i] = c;
    }
    fclose (fp);
  }
  printf ("curve=%s\nt=%d\n", curveInfo[curve].name, t);
  printDiagnostics (&bs, &st);
  return 0;
}
