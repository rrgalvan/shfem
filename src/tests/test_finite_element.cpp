#include <shfem.hpp>

// Boost test suite library, http://www.boost.org/doc/libs/release/libs/test/.
// Enable dynamic linking with -lboost_unit_test_framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE BasicTests
#include <boost/test/unit_test.hpp>

typedef shfem::real_t Real;
typedef shfem::dim2::Point Point;
typedef shfem::dim2::TriangleMesh Mesh;
typedef shfem::dim2::QuadRule<shfem::VerticesQR> VerticesQuadRule;

BOOST_AUTO_TEST_CASE(ElementAreaInOrder)
{
  // Define mesh
  Mesh mesh;
  mesh.read_file_msh("squared-mesh-2x2.msh");

  // Define finite element and attach to it the quadrature rule
  shfem::FiniteElement fe;

  // For every mesh cell...
  for(unsigned i=0; i<mesh.get_ncel(); ++i) {
    fe.reinit(mesh, i);  // Compute element-specific data for cell i
    // Compute area of current cell
    Real current_area = fe.area(), expected_area = 0.125;
    // Check that result differs from exact value no more than $10^{-20}$
    BOOST_CHECK_CLOSE(current_area, expected_area, 1.e-20);
  }
}

BOOST_AUTO_TEST_CASE(FiniteElementInOrder)
{
  // Define mesh
  Mesh mesh;
  mesh.read_file_msh("squared-mesh-2x2.msh");

  // Define quadrature rule
  VerticesQuadRule qr;

  // Define finite element and attach to it the quadrature rule
  shfem::FiniteElement fe;
  fe.set_quadrature_rule(qr);

  //,-----------------------------------------------
  //| Test affine transformation (for every element)
  //`-----------------------------------------------
  for (size_t icell=0; icell<mesh.get_ncel(); ++icell)
    {
      fe.reinit(mesh, icell); // Compute element-specific data

      // Store reference points, *in positive order*
      std::vector<Point> hatP = {Point(0,0), Point(1,0), Point(0,1)};
      // Check that they are respectively transformed into cell vertices
      for (size_t ipoint = 0; ipoint < 3; ++ipoint)
	{
	  Point P0 = fe.apply_affine_map(hatP[ipoint]);
	  Point P1 = fe.get_vertex(ipoint);
	  BOOST_CHECK_EQUAL(P0.x, P1.x);
	  BOOST_CHECK_EQUAL(P0.y, P1.y);
	}
    }

  //,-------------------------------------------------------
  //| Test inverse affine transformation (for every element)
  //`-------------------------------------------------------
  for (size_t icell=0; icell<mesh.get_ncel(); ++icell)
    {
      // Store reference points, *in positive order*
      std::vector<Point> hatP = {Point(0,0), Point(1,0), Point(0,1)};
      fe.reinit(mesh, icell); // Compute element-specific data
      // Check that pysical verties are transformed into reference vertices
      for (size_t ipoint = 0; ipoint < 3; ++ipoint)
	{
	  Point P0 = hatP[ipoint];
	  Point P1 = fe.apply_inv_affine_map(fe.get_vertex(ipoint));
	  BOOST_CHECK_EQUAL(P0.x, P1.x);
	  BOOST_CHECK_EQUAL(P0.y, P1.y);
	}
    }
}

BOOST_AUTO_TEST_CASE(QuadRuleInOrder)
{
  int nb_nodes = 3;

  // Define a quadrature rule
  VerticesQuadRule qr;
  BOOST_CHECK_EQUAL(qr.size(), nb_nodes);

  // Define a function on quadrature rule nodes
  shfem::QuadFunction f(nb_nodes);
  for(int i=0; i<nb_nodes; i++) f[i]=1.0;  // Set f={1,1,1}

  // Assure integral on reference triangle of previous function is correct
  Real integral = qr.integrate_on_ref_element(f);
  Real area_of_ref_triangle = 0.5; // Area of reference triangle
  BOOST_CHECK_CLOSE(integral, area_of_ref_triangle, 1.e-20);
}
