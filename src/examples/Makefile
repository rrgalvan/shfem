# gcc compile
GCC = g++
GCFLAGS =-O3 -Wall -std=c++11 #-fopenmp -fPIC -shared

SHFEM_SRCDIR=..
EIGEN_DIR = /usr/include/eigen3/
INC = -I$(SHFEM_SRCDIR) -I$(EIGEN_DIR)

TARGET = read_msh integrate equation_system speedup-test

all: $(TARGET)

%: %.cpp
	$(GCC) $(GCFLAGS) $(INC) $< -o $@

run: $(TARGET)
	@for program in $(TARGET); do \
		echo "========================================"; \
		echo "Running $$program..."; \
		echo "========================================"; \
		./$$program ; \
	done

clean:
	rm -rf *.a *.o $(TARGET)
