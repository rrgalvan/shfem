HPFEM_SRCDIR=..

ALL: test_read_msh test_fespace

#REMARK: 
# It is important to locate "-lboost_unit_test_framework" at end of line,
# see http://askubuntu.com/questions/104316/problem-with-executing-boost-library

test_read_msh: test_read_msh.cpp
	g++ -I$(HPFEM_SRCDIR) $< -o $@ -lboost_unit_test_framework 

test_fespace: test_fespace.cpp
	g++ -I$(HPFEM_SRCDIR) $< -o $@ -lboost_unit_test_framework 
