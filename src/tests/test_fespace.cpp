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
  shfem::CG_FESpace fe(m);
  fe.get_mesh();

  shfem::dim2::QuadRule<shfem::VerticesQR> Q;
}
