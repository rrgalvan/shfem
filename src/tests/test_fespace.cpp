#include <shfem.hpp>

// Boost test suite library.
// Enable dynamic linking with -lboost_unit_test_framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE BasicTests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(FESpaceInOrder)
{
  shfem::dim2::Mesh m;
  m.read_file_msh("squared-mesh-2x2.msh");
  shfem::CG_FESpace fe(m);
  fe.get_mesh();
}

BOOST_AUTO_TEST_CASE(QuadRuleInOrder)
{
  typedef shfem::real_t real_t;
  shfem::dim2::QuadRule<shfem::VerticesQR> qr;

  int nvertices = 3;
  shfem::QuadFunction f(nvertices);
  for(int i=0; i<nvertices; i++) f[i]=1.0;  // Set f={1,1,1}

  real_t integral = qr.integrate_on_ref_element(f);
  real_t area_of_ref_triangle = 0.5; // Area of reference triangle
  BOOST_CHECK_CLOSE(integral, area_of_ref_triangle, 1.e-20);
}
