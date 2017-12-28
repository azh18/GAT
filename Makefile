# Makefile example for compiling cuda and linking cuda to cpp:

SOURCELOC = 

UTILITYLOC =

NEWMOD =

PROGRAM = CUDAGTS

INCDIR= .
#
# Define the C compile flags
CCFLAGS = -g -m64 -I /usr/local/cuda/include -I ./header -O2 -std=c++11
CC = g++

# Define the Cuda compile flags
#
CUDAFLAGS= -O2 --gpu-architecture=compute_37 --gpu-code=compute_37 -m64 -use_fast_math -I ./header -std=c++11 -lineinfo --use-local-env -ccbin "g++" -cudart static --cl-version 2015

CUDACC= nvcc

# Define Cuda objects

#

CUDA = kernel.o

# Define the libraries

SYSLIBS= -lc
USRLIB  = -lcudart

# Define all object files

OBJECTS = \
	CUDAGTS_CPP.o\
	Cell.o\
	DateTime.o\
	FVTable.o\
	Grid.o\
	MBB.o\
	PreProcess.o\
	QueryResult.o\
	SamplePoint.o\
	Trajectory.o\
	SystemTest.o\
	FSG.o\
	STIG.o\
	MortonGrid.o



install: CUDAGTS


# Define Task Function Program


all: CUDAGTS


# Define what Modtools is


CUDAGTS: $(OBJECTS) $(CUDA)

	$(CUDACC) $(CUDAFLAGS) -o CUDAGTS -L/usr/local/cuda/lib64 -lcuda $(OBJECTS) $(CUDA) $(USRLIB) $(SYSLIBS)

# Modtools_Object codes


CUDAGTS_CPP.o: main.cpp

	$(CC) $(CCFLAGS) -c main.cpp -o CUDAGTS_CPP.o

.cpp.o:

	$(CC) $(CCFLAGS) -c $<


CUDAINCDIR= /usr/local/cuda/include 

kernel.o: kernel.cu

	$(CUDACC) $(CUDAFLAGS) --compile -c kernel.cu -o kernel.o

clean:
	rm -rf *.o
#  end
