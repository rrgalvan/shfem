# gcc compiler
GCC = g++
GCFLAGS =-O3 -Wall -std=c++11 #-fopenmp -fPIC -shared

SHFEM_SRCDIR = ..

BOOST_UNIT_TEST_LIB = -lboost_unit_test_framework

TARGET = test_read_msh test_finite_element

all: $(TARGET)

#REMARK:
# It is important to locate "-lboost_unit_test_framework" at end of line,
# see http://askubuntu.com/questions/104316/problem-with-executing-boost-library

test_read_msh: test_read_msh.cpp
	$(GCC) $(GCFLAGS) -I$(SHFEM_SRCDIR) $< -o $@ $(BOOST_UNIT_TEST_LIB)

test_finite_element: test_finite_element.cpp
	$(GCC) $(GCFLAGS) -I$(SHFEM_SRCDIR) $< -o $@ $(BOOST_UNIT_TEST_LIB)

run: $(TARGET)
	@for program in $(TARGET); do \
		echo "========================================"; \
		echo "Running $$program..."; \
		echo "========================================"; \
		./$$program ; \
	done

clean:
	rm -rf *.a *.o $(TARGET)
