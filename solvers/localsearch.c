#include <stdlib.h>

#include "Instance.h"
#include "Utilities.h"
#include "tour.h"
#include "neighbors.h"
#include "localsearch.h"

/*********************************************************************

Local search: keep making small changes that shorten the tour until none
is left.  Two kinds of change ("moves"):

2-OPT  Remove two edges and reconnect the tour the other way.  This is what
       untangles two crossing edges:

	   a -> b                      a    b
	          \/        becomes    |    |
	          /\                   |    |
	   d <- c                      c    d          (path b..c is reversed)

       removed (a,b) + (c,d), added (a,c) + (b,d).

OR-OPT Cut a short chain of 1-3 consecutive cities out of the tour, close
       the gap, and re-insert the chain (either way round) between two other
       neighbouring cities u and v:

	   a [s1 .. s2] b   ...   u v     becomes     a b  ...  u [s1 .. s2] v

Speed tricks (all standard, see Bentley 1992; Johnson & McGeoch 1997):

  * Neighbour lists: a move must connect a city to one of its nbK nearest
    neighbours (neighbors.c), so each search is O(nbK), not O(N).
  * Work queue ("don't-look bits"): only cities in the queue are searched.
    Initially every city is queued.  A city that yields no improving move
    leaves the queue; whenever a move changes the edges at a city, that city
    is queued again.  The search ends when the queue is empty, i.e. the tour
    is a local optimum: no single 2-opt (or Or-opt) move can improve it.
  * Reversal cap: a move that would reverse more than MAX_REVERSAL cities
    is skipped.  On an array tour each reversal costs time proportional to
    its length, and on million-city tours unlimited reversals make the
    search roughly O(N^2).  Tours up to 2*MAX_REVERSAL cities are unaffected
    (the shorter side is never longer than half the tour).

Moves are applied as one or more 2-opt exchanges (tour.c exchange2), which
reverse the shorter side of the tour.

 *********************************************************************/

#define D(a,b) distance ((a), (b))

#ifndef MAX_REVERSAL          /* can be overridden: make CFLAGS="-O2 -DMAX_REVERSAL=..." */
#define MAX_REVERSAL 100000   /* skip moves that reverse more cities than this */
#endif

/* ------------------------------------------------------------ work queue */

static int *queue;            /* circular FIFO of cities to (re)examine */
static int *queued;           /* queued[c] = 1 while c is in the queue  */
static int  qHead, qCount;

static void push (int c) {

  int tail;

  if (queued[c]) return;                 /* already waiting: don't add twice */
  tail = qHead + qCount;                 /* slot after the last element ...  */
  if (tail >= ncities) tail -= ncities;  /* ... wrapping round the array     */
  queue[tail] = c;
  queued[c] = 1;
  qCount++;
}

static int pop (void) {

  int c = queue[qHead];                  /* take the oldest city */
  if (++qHead == ncities) qHead = 0;
  qCount--;
  queued[c] = 0;
  return c;
}

/* ------------------------------------------------------------ 2-opt */

/* Look for the best 2-opt move that removes one of the two tour edges at
   city a and connects a to one of its neighbours c.  Apply it if it
   shortens the tour.  Returns 1 if a move was made. */

static int twoOptMove (int a) {

  int dir, b, c, d, k, g1, g, delta;
  int best = 0, ba = 0, bb = 0, bc = 0, bd = 0;
  int *nb = neighbors (a);

  for (dir = 0; dir < 2; ++dir) {        /* try both tour edges at a:        */
    b = dir ? pred(a) : succ(a);         /*   dir 0: a -> b,  dir 1: b -> a  */
    g1 = D(a,b);                         /* length we would gain by removing (a,b) */

    for (k = 0; k < nbK; ++k) {          /* candidate new partners, nearest first */
      c = nb[k];
      g = g1 - D(a,c);                   /* gain so far: remove (a,b), add (a,c)  */
      if (g <= 0) break;                 /* farther neighbours only get worse:    */
                                         /*   stop (standard gain criterion)      */
      d = dir ? pred(c) : succ(c);       /* c's tour edge pointing the same way   */
      if (c == b || d == a) continue;    /* edges touch: not a real 2-opt move    */

      delta = g + D(c,d) - D(b,d);       /* full saving: also remove (c,d),       */
                                         /*   add (b,d)                           */
      if (delta > best                   /* best move found so far, and not  */
	  && exchangeCost (a, b, c, d) <= MAX_REVERSAL) {  /* too slow: keep it */
	best = delta;
	ba = a; bb = b; bc = c; bd = d;
      }
    }
  }

  if (best <= 0) return 0;               /* no improving move at this city */

  exchange2 (ba, bb, bc, bd);            /* remove (a,b),(c,d); add (a,c),(b,d) */
  push (ba); push (bb);                  /* their edges changed: look again     */
  push (bc); push (bd);
  return 1;
}

/* ------------------------------------------------------------ Or-opt */

/* Look for the best Or-opt move that moves a chain of 1-3 cities starting
   at city s1 to a better place.  Returns 1 if a move was made.

   "Forward" below means along the direction dir: dir 0 follows succ(),
   dir 1 follows pred().  In that direction the tour reads

	   a  s1 .. s2  b  .....  u  v  .....        (v comes right after u)

   and the move turns it into either

	   a  b  .....  u  s1 .. s2  v  .....        "same way round", or
	   a  b  .....  u  s2 .. s1  v  .....        "reversed".              */

static int orOptMove (int s1) {

  int dir, len, s2, mid, a, b, g0, e, k, c, side, u, v, dSame, dRev;
  int best = 0, ba = 0, bs1 = 0, bs2 = 0, bb = 0, bu = 0, bv = 0;
  int bsame = 0;
  int *nb;

  for (dir = 0; dir < 2; ++dir) {                  /* chains going both ways */
    s2 = s1;
    mid = s1;
    for (len = 1; len <= 3; ++len) {               /* chain lengths 1, 2, 3 */
      if (len > 1) {                               /* extend chain by one city */
	mid = s2;                                  /*   (middle city if len 3) */
	s2 = dir ? pred(s2) : succ(s2);
      }
      a = dir ? succ(s1) : pred(s1);               /* city before the chain */
      b = dir ? pred(s2) : succ(s2);               /* city after the chain  */
      if (a == b || a == s2) break;                /* tour too small        */

      g0 = D(a,s1) + D(s2,b) - D(a,b);             /* saving from cutting the */
      if (g0 <= 0) continue;                       /*   chain out & closing   */
                                                   /*   the gap with (a,b)    */

      for (e = 0; e < 2; ++e) {                    /* chain end that gets the  */
	int end = e ? s2 : s1;                     /*   new neighbour: s1, s2  */
	if (e == 1 && len == 1) break;             /* one-city chain: same end */
	nb = neighbors (end);

	for (k = 0; k < nbK; ++k) {                /* candidate insertion spots */
	  c = nb[k];
	  if (D(end,c) >= g0) break;               /* new edge already eats the */
	                                           /*   whole saving: stop      */
	  for (side = 0; side < 2; ++side) {       /* insert after c or before c */
	    if (side == 0) { u = c;  v = dir ? pred(c) : succ(c); }
	    else           { v = c;  u = dir ? succ(c) : pred(c); }

	    if (u == s1 || u == s2 || u == mid ||  /* (u,v) must not touch */
		v == s1 || v == s2 || v == mid)    /*   the chain itself   */
	      continue;

	    /* total saving = cut-out saving + edge (u,v) removed
	       - the two edges that join the chain to u and v */
	    dSame = g0 + D(u,v) - D(u,s1) - D(s2,v);   /* u s1..s2 v */
	    dRev  = g0 + D(u,v) - D(u,s2) - D(s1,v);   /* u s2..s1 v */

	    if ((dSame > best || dRev > best)          /* best so far, and the */
		&& tourGap (s1, u) <= MAX_REVERSAL) {  /*   chain doesn't move */
	                                               /*   too far: keep it   */
	      best  = (dSame > dRev) ? dSame : dRev;
	      bsame = (dSame > dRev);
	      ba = a; bs1 = s1; bs2 = s2; bb = b; bu = u; bv = v;
	    }
	  }
	}
      }
    }
  }

  if (best <= 0) return 0;                         /* nothing improves */

  /* Apply the move as a sequence of 2-opt exchanges.  Each line removes
     two edges and adds two (see exchange2 in tour.c); after the first one
     or two lines the chain sits between u and v reversed (u s2..s1 v). */

  if (bv == ba) {                    /* u a [s1..s2] b  ->  u [s2..s1] a b       */
    exchange2 (bu, ba, bs2, bb);     /*   remove (u,a),(s2,b); add (u,s2),(a,b)  */
  }
  else if (bu == bb) {               /* a [s1..s2] b v  ->  a b [s2..s1] v       */
    exchange2 (ba, bs1, bb, bv);     /*   remove (a,s1),(b,v); add (a,b),(s1,v)  */
  }
  else {                             /* general case, two steps:                 */
    exchange2 (ba, bs1, bu, bv);     /*   a [s1..s2] b..u v -> a u..b [s2..s1] v */
    exchange2 (ba, bu, bb, bs2);     /*   a u..b [s2..s1] v -> a b..u [s2..s1] v */
  }
  if (bsame && bs1 != bs2)           /* wanted "same way round": flip the chain  */
    exchange2 (bu, bs2, bs1, bv);    /*   u [s2..s1] v -> u [s1..s2] v           */

  push (ba); push (bb);              /* every city whose edges changed */
  push (bs1); push (bs2);
  push (bu); push (bv);
  return 1;
}

/* ------------------------------------------------------------ driver */

void localSearch (int useOrOpt, long long *twoOptMoves, long long *orOptMoves) {

  int i, c;

  *twoOptMoves = 0;
  *orOptMoves = 0;
  if (ncities < 8) return;               /* tiny tours: moves would overlap */

  queue  = (int *) malloc (ncities * sizeof(int));
  queued = (int *) calloc (ncities + 1, sizeof(int));
  if (queue == 0 || queued == 0) fatal ("Out of memory in localSearch");
  qHead = 0;
  qCount = 0;

  for (i = 0; i < ncities; ++i) push (tour[i]);   /* start: every city queued */

  while (qCount > 0) {                   /* until no city has work left:     */
    c = pop();                           /*   take the next city             */
    if (twoOptMove (c)) {                /*   try 2-opt first                */
      (*twoOptMoves)++;
      push (c);                          /*   improved: check c again later  */
    }
    else if (useOrOpt && orOptMove (c)) {/*   else try Or-opt                */
      (*orOptMoves)++;
      push (c);
    }
    /* no move found: c stays out of the queue until a neighbour's
       move changes its edges and pushes it back in */
  }

  free (queue);
  free (queued);
}
