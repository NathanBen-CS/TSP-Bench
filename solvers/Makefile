# Makefile for your own TSP solvers.
#
# Every solver is one .c file with a main(), linked with the DIMACS reader
# (../codes/Instance.c) and helpers (../codes/Utilities.c).  The local
# search solvers also link the shared pieces in LIBOBJS.
#
# To add a solver "foo":  write foo.c, add "foo" to SOLVERS (and to LS_SOLVERS
# if it uses the shared library), run make.
#
# Keep CFLAGS identical for every solver: the benchmark rules require the
# same settings everywhere, and bench/machine_bench.py records this line.

CC = gcc
CFLAGS = -O2
DIMACS = ../codes

LS_SOLVERS = spacefill twoopt twoopt_or
SOLVERS = nn $(LS_SOLVERS)

LIBOBJS = tour.o neighbors.o hilbert.o greedystart.o localsearch.o
HEADERS = tour.h neighbors.h hilbert.h greedystart.h localsearch.h \
	$(DIMACS)/Instance.h $(DIMACS)/Utilities.h

all:	$(SOLVERS)

Instance.o: $(DIMACS)/Instance.c $(DIMACS)/Instance.h $(DIMACS)/Utilities.h
	$(CC) $(CFLAGS) -I$(DIMACS) -c $(DIMACS)/Instance.c -o Instance.o

Utilities.o: $(DIMACS)/Utilities.c $(DIMACS)/Utilities.h
	$(CC) $(CFLAGS) -I$(DIMACS) -c $(DIMACS)/Utilities.c -o Utilities.o

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -I$(DIMACS) -c $< -o $@

nn: nn.o Instance.o Utilities.o
	$(CC) $(CFLAGS) $^ -lm -o $@

$(LS_SOLVERS): %: %.o $(LIBOBJS) Instance.o Utilities.o
	$(CC) $(CFLAGS) $^ -lm -o $@

clean:
	rm -f *.o *.exe $(SOLVERS) *~
