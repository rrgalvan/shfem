#include <shfem.hpp>

// Boost test suite library.
// Enable dynamic linking with -lboost_unit_test_framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE BasicTests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(MeshReadInOrder)
{
  shfem::TriangleMesh m;
  m.read_file_msh("squared-mesh-2x2.msh");

  unsigned nb_triangles = m.get_ncel(), expected_nb_triangles = 8;
  BOOST_CHECK_EQUAL(nb_triangles, expected_nb_triangles);

  unsigned nb_vertices = m.get_nver(), expected_nb_vertices = 9;
  BOOST_CHECK_EQUAL(nb_vertices, expected_nb_vertices);
}
