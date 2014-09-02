#include <hpfem.hpp>

// Boost test suite library.
// Enable dynamic linking with -lboost_unit_test_framework
#define BOOST_TEST_DYN_LINK 
#define BOOST_TEST_MODULE BasicTests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(MeshReadInOrder)
{
  hpfem::Mesh m;
  m.read_file_msh("squared-mesh-2x2.msh");
  m.print();
  
  for(unsigned i=0; i<m.get_nt(); ++i) {
    hpfem::real_t
      current_area = m.area(i), 
      expected_area = 0.5;
    std::cout << "Area of triangle " << i << ": "
	      << current_area << std::endl; 
    
    // Check that floating point values differ no more 
    // than 0.0000000001 of their value.
    BOOST_CHECK_CLOSE(current_area, expected_area, 1.e-10);
  }
}
