# Makefile for your own TSP solvers.
#
# Every solver is one .c file here, linked with the DIMACS reader
# (../codes/Instance.c) and helpers (../codes/Utilities.c).
#
# To add a solver "foo":  write foo.c, add "foo" to SOLVERS below, run make.
#
# Keep CFLAGS identical for every solver: the benchmark rules require the
# same settings everywhere, and bench/machine_bench.py records this line.

CC = gcc
CFLAGS = -O2
DIMACS = ../codes

SOLVERS = nn

all:	$(SOLVERS)

Instance.o: $(DIMACS)/Instance.c $(DIMACS)/Instance.h $(DIMACS)/Utilities.h
	$(CC) $(CFLAGS) -I$(DIMACS) -c $(DIMACS)/Instance.c -o Instance.o

Utilities.o: $(DIMACS)/Utilities.c $(DIMACS)/Utilities.h
	$(CC) $(CFLAGS) -I$(DIMACS) -c $(DIMACS)/Utilities.c -o Utilities.o

%.o: %.c $(DIMACS)/Instance.h $(DIMACS)/Utilities.h
	$(CC) $(CFLAGS) -I$(DIMACS) -c $< -o $@

$(SOLVERS): %: %.o Instance.o Utilities.o
	$(CC) $(CFLAGS) $^ -lm -o $@

clean:
	rm -f *.o *.exe $(SOLVERS) *~
