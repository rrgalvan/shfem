// equation_system.cpp ---

// Copyright (C) 2015 Rafa Rodríguez Galván <rafaelDOTrodriguezATucaDOTes>

// Author: Rafa Rodríguez Galván

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include <chrono>  // std time handling (C++11)
#include <shfem.hpp>
#include <cmath>
#include <fstream> // To save the solution to a file

using namespace shfem;

#include <Eigen/Dense>
// using namespace Eigen;

Real rhs_function(Real x, Real y) {
  // Right hand side function
  return 4; // Exact solution: u = 1-x^2-y^2 in Omega=unit_circle
  // const double PI = 3.1415926535;
  // return 2*PI*PI*sin(PI*x)*sin(PI*y);
}

Real homog_dirichlet_function(Real x, Real y) {
  // Homogeneos Dirichlet condition
  return 0.0;
}

int main()
{
  // Try to read mesh contents from a .msh file (for information about the
  // structure of these files, see FreeFem++ documentation).
  TriangleMesh mesh;
  try {
    // mesh.read_file_msh("circle200.msh");
    mesh.read_file_msh("squared-mesh-2x2.msh"); }
  catch (...) {
    std::cerr << "Error reading mesh file" << std::endl;
    exit(1); }

  // Define quadrature rule with nodes located at vertices of the reference triangle
  VerticesQuadRule quad_rule;

  // Define a finite element space on current mesh
  P1_FE_Space fe_space(mesh, quad_rule);

  // Global finite element matrix
  int N = fe_space.get_ndofs();
  Eigen::MatrixXf A(N, N);

  // Global finite element rhs vector
  Eigen::VectorXf b(N);

  // Start chronometer
  auto start = std::chrono::high_resolution_clock::now();

  // For each cell, r:
  for (Index r=0; r<mesh.get_ncel(); ++r)
    {
      auto fe = fe_space.get_element(r); // Build a finite element on cell r

      // Get the basis functions of curent element
      const std::vector<FE_Function>& phi = fe.get_phi();
      // Get x-derivatives of all the basis functions of curent element
      const std::vector<FE_Function>& dx_phi = fe.get_dx_phi();
      // Get also y-derivatives of basis functions
      const std::vector<FE_Function>& dy_phi = fe.get_dy_phi();

      // // BORRAR
      // cout << endl;
      // for(int i=0; i<2; i++) cout << "," << dx_phi[0][i];
      // cout << endl;
      // for(int i=0; i<2; i++) cout << "," << dy_phi[0][i];

      // cout << endl;
      // for(int i=0; i<2; i++) cout << "," << dx_phi[1][i];
      // cout << endl;
      // for(int i=0; i<2; i++) cout << "," << dy_phi[1][i];

      // cout << endl;
      // for(int i=0; i<2; i++) cout << "," << dx_phi[2][i];
      // cout << endl;
      // for(int i=0; i<2; i++) cout << "," << dy_phi[2][i];

      // cout << endl;

      // Local stiffness matrix
      Index ndofs = fe.get_ndofs();
      Eigen::MatrixXf A_r(ndofs,ndofs);

      // Local rhs vector
      Eigen::VectorXf b_r(ndofs);

      // For each degree of freedom, ii
      for (Index i = 0; i < ndofs; ++i)
	{
	  // For each degree of freedom, j
	  for (Index j = 0; j < ndofs; ++j)
	    {
	      // Compute integral of product of x-derivatives
	      Real Ax_ij = fe.integrate(dx_phi[i], dx_phi[j]);

	      // Compute integral of product of y-derivatives
	      Real Ay_ij = fe.integrate(dy_phi[i], dy_phi[j]);

	      // Store the result in the local matrix
	      A_r(i,j) = Ax_ij + Ay_ij;
	    }

	  // Store also the i-th element in the rhs vector
	  auto qrule = fe.get_quadrature_rule();
	  FE_Function f(qrule.size()); // FE_Function on quadrature nodes
	  for(Index r=0; r<qrule.size(); r++) {
	    const FiniteElement::POINT& P = qrule.nodes[r];
	    f[i] = rhs_function(P.x, P.y);
	  }
	  b_r(i) = fe.integrate(f, phi[i]);
	}

      // Add local matrix A_r into global matrix A
      fe_space.add_local_matrix(A_r, A, r);

      // Add local vector b_r into global rhs vector b
      fe_space.add_local_vector(b_r, b, r);

    }

  // Print matrix and vector assembiling time
  auto end_system_time = std::chrono::high_resolution_clock::now();
  int elapsed_time =
    std::chrono::duration_cast<std::chrono::milliseconds>( end_system_time - start ).count();
  std::cout << "Equations system assembling time: " << elapsed_time << " miliseconds" << std::endl;

  // Apply Diriclet Conditions on matrix and vector
  DirichletConditions dirichlet;
  dirichlet[1] = homog_dirichlet_function;  // Condition on boundary 1
  dirichlet[2] = homog_dirichlet_function;  // Condition on boundary 1
  dirichlet[3] = homog_dirichlet_function;  // Condition on boundary 1
  dirichlet[4] = homog_dirichlet_function;  // Condition on boundary 1
  fe_space.apply_dirichlet_conditions(dirichlet, A, b);

  // Print Dirichlet conditions assembling time
  auto end_dirichlet_time = std::chrono::high_resolution_clock::now();
  elapsed_time =
    std::chrono::duration_cast<std::chrono::milliseconds>( end_dirichlet_time -
  							   end_system_time ).count();
  std::cout << "Dirichlet conditions assembling time: " << elapsed_time << " miliseconds" << std::endl;

  // Print equations system
  std::cout << "Resulting matrix:" << std::endl << A << std::endl;
  std::cout << "Resulting vector:" << std::endl << b << std::endl;
  Eigen::VectorXf u =  A.partialPivLu().solve(b);
  std::cout << "Solution: " << u << std::endl;

  // Save solution to file (so that it can be read by external tools
  // like FreeFem++)
  std::ofstream myfile;
  std::string sol_filename = "solution.txt";
  myfile.open (sol_filename);
  std::cout << "Writing solution to file: " << sol_filename << endl;
  myfile << u.size() << endl << u;
  myfile.close();
  return 0;
}
