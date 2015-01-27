#include <shfem.hpp>

// Boost test suite library, http://www.boost.org/doc/libs/release/libs/test/.
// Enable dynamic linking with -lboost_unit_test_framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE BasicTests
#include <boost/test/unit_test.hpp>

typedef shfem::Real Real;
typedef shfem::Index Index;
typedef shfem::dim2::Point Point;
typedef shfem::dim2::TriangleMesh Mesh;
typedef shfem::dim2::QuadRule<shfem::VerticesQR> VerticesQuadRule;

BOOST_AUTO_TEST_CASE(ElementAreaTest)
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

inline Index delta_kronecker(Index i, Index j) {
  return (i==j) ? 1 : 0;
}

BOOST_AUTO_TEST_CASE(BasisFunctionsTest)
{
  // Define mesh
  Mesh mesh;
  mesh.read_file_msh("squared-mesh-2x2.msh");

  // Define finite element and attach to it the quadrature rule
  shfem::FiniteElement fe;
  // Define quadrature rule.
  VerticesQuadRule qr;
  fe.set_quadrature_rule(qr);

  for (Index icell=0; icell<mesh.get_ncel(); ++icell)
    {
      fe.reinit(mesh, icell); // Compute element-specific data
      // For each degree of freedom, i:
      for (Index i = 0; i<fe.get_ndofs(); ++i)
	{
	  // Compute i-th basis function, evaluated at quadrature nodes
	  const shfem::FE_Function& phi_i = fe.get_basis_function(i);
	  // For each quadrature point, j:
	  for(Index j=0; j<qr.size(); ++j)
	    // For this quadrture rule, nodes are equal to dof
	    // (vertices), therefore phi_i(node_j) = delta_kronecker(i,j)
	    BOOST_CHECK_EQUAL(phi_i[j], delta_kronecker(i,j));
	}
    }
}

BOOST_AUTO_TEST_CASE(AffineMapTest)
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
  for (Index icell=0; icell<mesh.get_ncel(); ++icell)
    {
      fe.reinit(mesh, icell); // Compute element-specific data

      // Store reference points, *in positive order*
      std::vector<Point> hatP = {Point(0,0), Point(1,0), Point(0,1)};
      // Check that they are respectively transformed into cell vertices
      for (Index ipoint = 0; ipoint < 3; ++ipoint)
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
  for (Index icell=0; icell<mesh.get_ncel(); ++icell)
    {
      // Store reference points, *in positive order*
      std::vector<Point> hatP = {Point(0,0), Point(1,0), Point(0,1)};
      fe.reinit(mesh, icell); // Compute element-specific data
      // Check that pysical verties are transformed into reference vertices
      for (Index ipoint = 0; ipoint < 3; ++ipoint)
	{
	  Point P0 = hatP[ipoint];
	  Point P1 = fe.apply_inv_affine_map(fe.get_vertex(ipoint));
	  BOOST_CHECK_EQUAL(P0.x, P1.x);
	  BOOST_CHECK_EQUAL(P0.y, P1.y);
	}
    }
}

BOOST_AUTO_TEST_CASE(QuadRuleTest)
{
  int nb_nodes = 3;

  // Define a quadrature rule
  VerticesQuadRule qr;
  BOOST_CHECK_EQUAL(qr.size(), nb_nodes);

  // Define a function on quadrature rule nodes
  shfem::FE_Function f(nb_nodes);
  for(int i=0; i<nb_nodes; i++) f[i]=1.0;  // Set f={1,1,1}

  // Assure integral on reference triangle of previous function is correct
  Real integral = qr.integrate_on_ref_element(f);
  Real area_of_ref_triangle = 0.5; // Area of reference triangle
  BOOST_CHECK_CLOSE(integral, area_of_ref_triangle, 1.e-20);
}

BOOST_AUTO_TEST_CASE(IntegrateTest)
{
  // Define mesh on (0,1)x(0,1) with 2 subintervals on (0,1)
  Mesh mesh;
  mesh.read_file_msh("squared-mesh-2x2.msh");

  // Define quadrature rule
  VerticesQuadRule qr;

  // Define finite element and attach to it the quadrature rule
  shfem::FiniteElement fe;
  fe.set_quadrature_rule(qr);

  for (Index icell=0; icell<mesh.get_ncel(); ++icell)
    {
      fe.reinit(mesh, icell); // Compute element-specific data
      // For each degree of freedom, i:
      for (Index i = 0; i<fe.get_ndofs(); ++i)
	{
	  // Constant function, f=1, on quadrature nodes
	  const shfem::FE_Function f = {1., 1., 1.};
	  shfem::Real integral_f_f = fe.integrate(f,f); // integral of f*f
	  BOOST_CHECK_EQUAL(integral_f_f, fe.area());
	  BOOST_CHECK_EQUAL(integral_f_f, 0.5*0.5/2.); // Area of triangle: b*h/2

	  // Basis function on dof i
	  const shfem::FE_Function& phi_i = fe.get_basis_function(i);
	  // If K = Triang(P0,P1,P2), P0=(0,0), P1=(0,0.5), P2=(0.5,0.5),
	  // then $\phi_0(x,y) = 1-2x$ and $\int_K 1-2*x dy dx =
	  // = \int_0^0.5 \int_0^y 1-2x dy dx = 1/24$
	  shfem::Real integral_f_phi_i = fe.integrate(f, phi_i);
	  BOOST_CHECK_EQUAL(integral_f_phi_i, 1./24);

	  // $\phi_0^2(x,y) = (1-2x)^2$ and $\int_K (1-2*x) dy dx =
	  // = \int_0^0.5 \int_0^y (1-2x) dy dx = 1/48$
	  // But this quadrature rule is not exact for order 2 polynomials!!
	  shfem::Real integral_phi_i_phi_i = fe.integrate(phi_i, phi_i);
	  BOOST_CHECK_NE(integral_phi_i_phi_i, 1./48);
	}
    }
}
