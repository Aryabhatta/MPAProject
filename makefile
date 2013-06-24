all: project clean

# which compiler
CC = g++

# Options while compiling
CFLAGS = -g 

# Options for CFITSIO
CFITSIO = -L. -lcfitsio -lm

# Options for FFTW
FFTW =  -lfftw3 -lfftw3_omp

project: getchi.o readGridfile.o ecorr.o cont_rscl.o get_rv.o lfp.o read_lf_grid.o modules.o idlFuncn.o
	$(CC) -o project getchi.o readGridfile.o ecorr.o cont_rscl.o get_rv.o lfp.o read_lf_grid.o modules.o idlFuncn.o $(CFITSIO) $(FFTW)

getchi.o: getchi.cpp
	$(CC) $(CFLAGS) -c getchi.cpp $(CFITSIO)
	
readGridfile.o: readGridfile.cpp readGridfile.hpp
	$(CC) $(CFLAGS) -c readGridfile.cpp
	
ecorr.o: ecorr.cpp ecorr.hpp
	$(CC) $(CFLAGS) -c ecorr.cpp
	
cont_rscl.o: cont_rscl.cpp cont_rscl.hpp
	$(CC) $(CFLAGS) -c cont_rscl.cpp 
	
get_rv.o: get_rv.cpp get_rv.hpp
	$(CC) $(CFLAGS) -c get_rv.cpp

lfp.o: lfp.cpp lfp.hpp
	$(CC) $(CFLAGS) -c lfp.cpp

read_lf_grid.o: read_lf_grid.cpp read_lf_grid.hpp
	$(CC) $(CFLAGS) -c read_lf_grid.cpp

modules.o: modules.cpp modules.hpp
	$(CC) $(CFLAGS) -c modules.cpp

idlFuncn.o: idlFuncn.cpp idlFuncn.hpp
	$(CC) $(CFLAGS) -c idlFuncn.cpp

# Cleaning the intermediates
clean:
#	rm -f *.o
	rm -f *~
	rm -f *.gch

# Cleaning all  the intermediates
cleanAll:
	rm -f *.o
	rm -f *~
	rm -f *.gch
