# gcc compile
GCC = g++
GCFLAGS =-O3 -Wall -std=c++11 #-fopenmp -fPIC -shared

SHFEM_SRCDIR=..

TARGET = read_msh

all: $(TARGET)

read_msh: read_msh.cpp
	$(GCC) $(GCFLAGS) -I$(SHFEM_SRCDIR) $< -o $@

clean: 
	rm -rf *.a *.o $(TARGET)