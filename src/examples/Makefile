# gcc compile
GCC = g++
GCFLAGS =-O3 -Wall -std=c++11 #-fopenmp -fPIC -shared

SHFEM_SRCDIR=..

TARGET = read_msh integrate

all: $(TARGET)

%: %.cpp
	$(GCC) $(GCFLAGS) -I$(SHFEM_SRCDIR) $< -o $@

run: $(TARGET)
	@for program in $(TARGET); do \
		echo "========================================"; \
		echo "Running $$program..."; \
		echo "========================================"; \
		./$$program ; \
	done

clean:
	rm -rf *.a *.o $(TARGET)
