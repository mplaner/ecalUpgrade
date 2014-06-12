CXX=g++
CC=gcc
CFLAGS=-O2 -Wall -Wno-deprecated-declarations
LDFLAGS=-Llib 
INS=-Iinclude 

CFLAGS += `root-config --cflags`
LIBS   += `root-config --glibs`

OBJS=analysis.o init.o root.o utils.o postInitialResultsAnalysis.o

.PHONY: clean all main test

all: $(OBJS)

main: main.o
	$(CXX) -o main main.o $(OBJS) $(LIBS)

test: all main
	./main

clean:
	@rm -f *.o *.exe core* main

##################### Rules #####################
.cc.o:
	$(CXX) $(CFLAGS) $(INS) -c $<

.cpp.o:
	$(CXX) $(CFLAGS) $(INS) -c $<
