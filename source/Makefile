#OBJS = keys.o readkey.o genutil.o mt19937.o inputparams.o networkutil.o  sort2int.o meshutil.o meshdynamics.o flowutil.o contractutil.o gmresutil.o dgmres.o zgmres.o

OBJS = keys.o readkey.o genutil.o mt19937.o inputparams.o networkutil.o  sort2int.o meshutil.o dynsysutil.o fvmdynamics.o nodecontractutil.o flowutil.o resvutil.o

FC = gfortran
FLAGS = -O3 
DEBUGFLAGS = -g -fbacktrace -ffpe-trap=zero,overflow,invalid -fcheck=all -fbounds-check
OSGFLAGS = -static
LIBS =  -llapack -lblas  -lgfortran -lm

# executable name
EXE = ../netmeshdynamicsFVM.exe
OSGEXE=../netmeshdynamicsFVM_osg.exe
TESTEXE = ../testrun.exe

.SUFFIXES: .o .f90 

# general suffix rules
%.o: %.f90 
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

%.o: tools/%.f*
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

default: $(OBJS) main.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) main.o -o $(EXE) $(LIBS)

debug: $(OBJS) main.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) main.o -o $(EXE) $(LIBS)

osg: $(OBJS) main.o
	$(FC) $(FLAGS) $(OSGFLAGS) $(OBJS) main.o -o $(OSGEXE) $(LIBS) 

test: $(OBJS) testmain.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) testmain.o -o $(TESTEXE) $(LIBS)

basic: $(OBJS) test.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) test.o -o $(TESTEXE) $(LIBS)

clean: 
	rm *.o *.mod

veryclean: clean
	rm *~ $(EXE)


# ------------
# individual file dependencies
# ------------
readkey.o: inputparams.o keys.o genutil.o mt19937.o
genutil.o: mt19937.o
networkutil.o: keys.o genutil.o
meshutil.o: networkutil.o resvutil.o
dynsysutil.o: meshutil.o keys.o
fvmdynamics.o: dynsysutil.o networkutil.o nodecontractutil.o keys.o
nodecontractutil.o: meshutil.o dynsysutil.o networkutil.o flowutil.o
