SYSTEM     = x86-64_osx
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /Applications/Ilog/CPLEX_Studio1262/cplex
CONCERTDIR    = /Applications/Ilog/CPLEX_Studio1262/concert

# ---------------------------------------------------------------------
# Compiler selection
# ---------------------------------------------------------------------

CCC = g++ -O0
CC  = gcc -O0
JAVAC = javac

# ---------------------------------------------------------------------
# Compiler options
# ---------------------------------------------------------------------

CCOPT = -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD
COPT  = -m64 -fPIC
JOPT  = -classpath $(CPLEXDIR)/lib/cplex.jar -O

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXJARDIR   = $(CPLEXDIR)/lib/cplex.jar
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -m64 -lm -lpthread -framework CoreFoundation -framework IOKit
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -m64 -lm -lpthread -lclapack -lcblas -framework CoreFoundation -framework IOKit
JAVA      = java -Djava.library.path=$(CPLEXDIR)/bin/x86_darwin8_gcc4.0 -classpath $(CPLEXJARDIR):

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CFLAGS  = $(COPT)  -I$(CPLEXINCDIR) 
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)
JCFLAGS = $(JOPT)

LPCC: main.o LPCCsolver.o CPLEXcontrol.o cutsgenerator.o rngs.o rvgs.o utilities.o
	$(CC) $(CFLAGS) main.o LPCCsolver.o CPLEXcontrol.o cutsgenerator.o rngs.o rvgs.o utilities.o -o LPCC $(CLNFLAGS)

main.o: main.c LPCC.h LPCCsolver.h rngs.h rvgs.h utilities.h
	$(CC) -c $(CFLAGS) main.c -o main.o
LPCCsolver.o: LPCCsolver.c LPCCsolver.h
	$(CC) -c $(CFLAGS) LPCCsolver.c -o LPCCsolver.o
CPLEXcontrol.o: CPLEXcontrol.c CPLEXcontrol.h
	$(CC) -c $(CFLAGS) CPLEXcontrol.c -o CPLEXcontrol.o
cutsgenerator.o: cutsgenerator.c cutsgenerator.h
	$(CC) -c $(CFLAGS) cutsgenerator.c -o cutsgenerator.o
rngs.o: rngs.c rngs.h
	$(CC) -c $(CFLAGS) rngs.c -o rngs.o
rvgs.o: rvgs.c rvgs.h
	$(CC) -c $(CFLAGS) rvgs.c -o rvgs.o
utilities.o: utilities.c utilities.h
	$(CC) -c $(CFLAGS) utilities.c -o utilities.o

