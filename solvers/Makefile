# Makefile for the TSP solvers, 2D and 3D.
#
#	make        build everything, both dimensions
#	make 2d     only the 2D programs        make 3d     only the 3D programs
#
# How the dimension is chosen: everything in common/ only uses distance()
# and the neighbour lists, so it is compiled ONCE (build/common) and linked
# into both.  The parts that need coordinates share a header in common/
# (neighbors.h, curves.h) and have one implementation per dimension:
#
#	            2D (build/d2)            3D (build/d3)
#	reader      ../codes/Instance.c      d3/instance3d.c
#	neighbours  d2/neighbors2d.c         d3/neighbors3d.c
#	curves      d2/curves2d.c, hilbert.c d3/curves3d.c
#
# Programs go to ../bin with a 2d / 3d suffix.  A main in common/ (nn, twoopt,
# twoopt_or, boxdiag, hkbound) becomes one program per dimension.
#
# Keep CFLAGS identical for every solver: the benchmark rules require the
# same settings everywhere, and bench/machine_bench.py records this line.

CC = gcc
CFLAGS = -O2
DIMACS = ../codes
BIN = ../bin
BC = ../build/common
B2 = ../build/d2
B3 = ../build/d3

SHARED = $(BC)/tour.o $(BC)/localsearch.o $(BC)/greedystart.o $(BC)/boxes.o \
	$(BC)/merges.o $(BC)/diag.o $(BC)/Utilities.o
LIB2 = $(SHARED) $(B2)/Instance.o $(B2)/neighbors2d.o $(B2)/curves2d.o $(B2)/hilbert.o \
	$(B2)/rings2d.o
LIB3 = $(SHARED) $(B3)/instance3d.o $(B3)/neighbors3d.o $(B3)/curves3d.o $(B3)/kdboxes3d.o

COMMON_MAINS = nn twoopt twoopt_or boxdiag hkbound
PROGS2 = $(COMMON_MAINS:%=$(BIN)/%2d) $(BIN)/spacefill2d $(BIN)/boxtsp2d $(BIN)/boxtest2d \
	$(BIN)/ringtest2d
PROGS3 = $(COMMON_MAINS:%=$(BIN)/%3d) $(BIN)/boxtsp3d $(BIN)/boxtest3d $(BIN)/unit3d \
	$(BIN)/length3d $(BIN)/portgen3d $(BIN)/portcgen3d

HEADERS = $(wildcard common/*.h) $(wildcard d2/*.h) $(wildcard d3/*.h) \
	$(DIMACS)/Instance.h $(DIMACS)/Utilities.h

all:	2d 3d
2d:	dirs $(PROGS2)
3d:	dirs $(PROGS3)
dirs:
	@mkdir -p $(BIN) $(BC) $(B2) $(B3)

.SECONDARY:

# ---- shared objects: compiled once, linked into both dimensions

$(BC)/%.o: common/%.c $(HEADERS)
	$(CC) $(CFLAGS) -I$(DIMACS) -Icommon -c $< -o $@

$(BC)/%.o: $(DIMACS)/%.c $(HEADERS)
	$(CC) $(CFLAGS) -I$(DIMACS) -c $< -o $@

# ---- 2D

$(B2)/Instance.o: $(DIMACS)/Instance.c $(HEADERS)
	$(CC) $(CFLAGS) -I$(DIMACS) -c $< -o $@

$(B2)/%.o: d2/%.c $(HEADERS)
	$(CC) $(CFLAGS) -I$(DIMACS) -Icommon -Id2 -c $< -o $@

$(B2)/%.o: ../tests/%.c $(HEADERS)
	$(CC) $(CFLAGS) -I$(DIMACS) -Icommon -Id2 -c $< -o $@

$(COMMON_MAINS:%=$(BIN)/%2d): $(BIN)/%2d: $(BC)/%.o $(LIB2)
	$(CC) $(CFLAGS) $^ -lm -o $@

$(BIN)/spacefill2d $(BIN)/boxtsp2d: $(BIN)/%2d: $(B2)/%.o $(LIB2)
	$(CC) $(CFLAGS) $^ -lm -o $@

# tests for the 2D pieces: run ../bin/boxtest2d (exit code 0 = all pass)
$(BIN)/boxtest2d: $(B2)/boxtest2d.o $(LIB2)
	$(CC) $(CFLAGS) $^ -lm -o $@

# tests for the hull rings: run ../bin/ringtest2d (exit code 0 = all pass)
$(BIN)/ringtest2d: $(B2)/ringtest2d.o $(LIB2)
	$(CC) $(CFLAGS) $^ -lm -o $@

# ---- 3D

$(B3)/%.o: d3/%.c $(HEADERS)
	$(CC) $(CFLAGS) -I$(DIMACS) -Icommon -Id3 -c $< -o $@

$(B3)/%.o: ../tests/%.c $(HEADERS)
	$(CC) $(CFLAGS) -I$(DIMACS) -Icommon -Id3 -c $< -o $@

$(COMMON_MAINS:%=$(BIN)/%3d): $(BIN)/%3d: $(BC)/%.o $(LIB3)
	$(CC) $(CFLAGS) $^ -lm -o $@

$(BIN)/boxtsp3d $(BIN)/boxtest3d $(BIN)/unit3d: $(BIN)/%: $(B3)/%.o $(LIB3)
	$(CC) $(CFLAGS) $^ -lm -o $@

# the DIMACS tour checker with the 3D reader, and the 3D generators
$(BIN)/length3d: $(BC)/length.o $(B3)/instance3d.o $(BC)/Utilities.o
	$(CC) $(CFLAGS) $^ -lm -o $@

$(BIN)/portgen3d $(BIN)/portcgen3d: $(BIN)/%: $(DIMACS)/%.c
	$(CC) $(CFLAGS) $< -lm -o $@

clean:
	rm -f $(BC)/*.o $(B2)/*.o $(B3)/*.o $(BIN)/*.exe

.PHONY: all 2d 3d dirs clean
