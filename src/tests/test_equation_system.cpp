#include <Eigen/Dense>
#include <shfem.hpp>
using namespace Eigen;
using namespace shfem;

// Boost test suite library.
// Enable dynamic linking with -lboost_unit_test_framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE BasicTests
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(StiffnessMatrixTest)
{

  // Try to read mesh contents from a .msh file (for information about the
  // structure of these files, see FreeFem++ documentation).
  TriangleMesh mesh;
  try {
    mesh.read_file_msh("squared-mesh-2x2.msh"); }
  catch (...) {
    std::cerr << "Error reading mesh file" << std::endl;
    exit(1); }

  // Define quadrature rule with nodes located at vertices of the reference triangle
  VerticesQuadRule quad_rule;

  // Define a finite element space on current mesh
  P1_FE_Space fe_space(mesh, quad_rule);

  // Define global matrix
  int N = fe_space.get_ndofs();
  BOOST_CHECK(N == 9);  // Number of DOFs in current mesh
  MatrixXf K(N, N);

  // For each cell, r:
  for (Index r=0; r<mesh.get_ncel(); ++r)
    {
      const FiniteElement & fe = fe_space.get_element(r);

      // Get x-derivatives of all the basis functions of curent element
      // (a FE_Function is a vector wich stores values at quadrature points)
      const std::vector<FE_Function>& dx_phi = fe.get_dx_phi();
      // Get also y-derivatives of basis functions
      const std::vector<FE_Function>& dy_phi = fe.get_dy_phi();

      Index ndofs = fe.get_ndofs();
      // Local stiffness matrix
      MatrixXf K_r(ndofs,ndofs);
      MatrixXf target_matrix(ndofs, ndofs);
      target_matrix << 1./4, -1./24, -1./8,
	-1./24, 1./24, 0,
	-1./8,  0,  1./8;

      // For each degree of freedom, i:
      for (Index i = 0; i < ndofs; ++i)
	{
	  // For each degree of freedom, j:
	  for (Index j = 0; j < ndofs; ++j)
	    {
	      // Compute integral on element r of gradient(phi_i)*gradient(phi_j),
	      Real Kx_ij = fe.integrate(dx_phi[i], dx_phi[j]);
	      Real Ky_ij = fe.integrate(dy_phi[i], dy_phi[j]);
	      K_r(i,j) = Kx_ij + Ky_ij;
	    }
	}
      // Check that all matrices match the target matrix
      BOOST_CHECK(K_r == target_matrix);

      // Add local matrix K_r into global matrix K
      fe_space.assemble_matrix(K_r, K);
    }
}
