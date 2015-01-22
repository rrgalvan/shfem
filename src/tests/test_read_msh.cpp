#include <shfem.hpp>

// Boost test suite library.
// Enable dynamic linking with -lboost_unit_test_framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE BasicTests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(MeshReadInOrder)
{
  shfem::dim2::Mesh m;
  m.read_file_msh("squared-mesh-2x2.msh");

  unsigned nb_triangles = m.get_nelt(), expected_nb_triangles = 8;
  BOOST_CHECK_EQUAL(nb_triangles, expected_nb_triangles);

  unsigned nb_vertices = m.get_nver(), expected_nb_vertices = 9;
  BOOST_CHECK_EQUAL(nb_vertices, expected_nb_vertices);

  for(unsigned i=0; i<m.get_nelt(); ++i) {
    shfem::real_t current_area = m.area(i), expected_area = 0.125;
    // Check that the two floating point values differ no more
    // than $10^{-20}$ of their value.
    BOOST_CHECK_CLOSE(current_area, expected_area, 1.e-20);
  }
}
