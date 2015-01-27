// integrate.cpp ---

// Copyright (C) 2014 Rafa Rodríguez Galván <rafaelDOTrodriguezATucaDOTes>

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

#include <shfem.hpp>

// Some convenient definitions (TODO: avoid necessity of these definitions)
typedef shfem::Real Real;
typedef shfem::Index Index;
typedef shfem::FE_Function FE_Function;
typedef shfem::TriangleMesh TriangleMesh;
typedef shfem::QuadRule<shfem::VerticesQR> VerticesQuadRule;

int main()
{
  // Try to read mesh contents from a .msh file (for information about the
  // structure of these files, see FreeFem++ documentation).
  TriangleMesh mesh;
  try
    {
      mesh.read_file_msh("squared-mesh-2x2.msh");
    }
  catch (...)
    {
      std::cerr << "Error reading mesh file" << std::endl;
      exit(1);
    }

  // Define generic finite element
  shfem::FiniteElement fe;

  // Use quadrature rule with nodes located at vertices of the reference triangle
  VerticesQuadRule qr;
  fe.set_quadrature_rule(qr);

  // For each cell, r:
  for (Index r=0; r<mesh.get_ncel(); ++r)
    {
      // Compute element-specific data for current cell, including values of basis
      // functions (and of its derivatives) on quadrature nodes
      fe.reinit(mesh, r);

      Index ndofs = fe.get_ndofs();
      // For each degree of freedom, i:
      for (Index i = 0; i<ndofs; ++i)
	{
	  std::cout << std::endl;

	  // Get x-derivative of i-th basis function (FE_Function is a vector
	  // wich stores values on each quadrature point)
	  const FE_Function& dx_phi_i = fe.get_dx_basis_function(i);
	  // Get also y-derivatives on each quadrature point
	  const FE_Function& dy_phi_i = fe.get_dy_basis_function(i);

	  // For each degree of freedom, j:
	  for (Index j = 0; j<ndofs; ++j)
	    {
	      std::cout << std::endl;
	      // Compute integral of product of x-derivatives
	      const FE_Function& dx_phi_j = fe.get_dx_basis_function(j);
	      Real Kx_ij = fe.integrate(dx_phi_i, dx_phi_j);
	      std::cout << "Kx[" << i << "][" << j << "]=" << Kx_ij << std::endl;

	      // Compute integral of product of y-derivatives
	      const FE_Function& dy_phi_j = fe.get_dy_basis_function(j);
	      Real Ky_ij = fe.integrate(dy_phi_i, dy_phi_j);
	      std::cout << "Ky[" << i << "][" << j << "]=" << Ky_ij << std::endl;

	      // Compute integral on element r of gradient(phi_i)*gradient(phi_j),
	      // i.e. dx(phi_i)*dx(phi_j) + dy(phi_i)*dy(phi_j)
	      Real K_ij = Kx_ij + Ky_ij;
	      std::cout << "K[" << i << "][" << j << "]=" << K_ij << std::endl;
	    }
	}
    }
}
