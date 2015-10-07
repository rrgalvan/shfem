# Definition of the compiler
FC=gfortran
OPTIONS=-Wall -Wpedantic -fbounds-check -fbacktrace -finit-real=nan
INCLUDE=-I../src
COMPILE=$(FC) $(INCLUDE) $(OPTIONS) -c
LINK=$(FC)

# Localization of HOFFElib
HOFFELIB_SRC=../src
HOFFELIB=$(HOFFELIB_SRC)/libhoffe.a

# Definition of object and source files

SOURCES:= $(wildcard *.f90)
OBJECTS:=$(SOURCES:.f90=.o)
EXEC:=$(SOURCES:.f90=)

# Compiling rules

#all: mesh1d_example
all: $(EXEC)

%.o: %.f90
	$(COMPILE) $<
%: %.o
	$(LINK) $<

%: %.o
	$(LINK) $(OBJECTS) $(HOFFELIB) -o $@

clean:
	rm -f  $(EXEC) *.o *.mod *~

.PHONY: all clean
