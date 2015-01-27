# gcc compiler
GCC = g++
GCFLAGS =-O3 -Wall -std=c++11 #-fopenmp -fPIC -shared

SHFEM_SRCDIR = ..
EIGEN_DIR = /usr/include/eigen3/
INC = -I$(SHFEM_SRCDIR) -I$(EIGEN_DIR)

BOOST_UNIT_TEST_LIB = -lboost_unit_test_framework

TARGET = test_read_msh test_finite_element test_equation_system

all: $(TARGET)

#REMARK:
# It is important to locate "-lboost_unit_test_framework" at end of line,
# see http://askubuntu.com/questions/104316/problem-with-executing-boost-library
%: %.cpp
	$(GCC) $(GCFLAGS) $(INC) $< -o $@ $(BOOST_UNIT_TEST_LIB)

run: $(TARGET)
	@for program in $(TARGET); do \
		echo "========================================"; \
		echo "Running $$program..."; \
		echo "========================================"; \
		./$$program ; \
	done

clean:
	rm -rf *.a *.o $(TARGET)
