// integrate.cpp ---

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
using namespace shfem;

int main()
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

  // Start chronometer
  auto start = std::chrono::high_resolution_clock::now();

  // Repeat several times
  int Nrepeat = 5000;
  for (int num=0; num<Nrepeat; ++num)
    {
      // For each cell, r:
      for (Index r=0; r<mesh.get_ncel(); ++r)
	{
	  // In C++11, one can use "auto" instead of "FiniteElement"

	  auto fe = fe_space.get_element(r); // Build a finite element on cell r

	  // Get x-derivatives of all the shape functions on curent element.
	  //
	  // The get_dx_phi() member function returns (a const reference
	  // to) a vector of FE_Function objects (one FE_Function for each
	  // dof). The i-th FE_Function represents (the evaluation on
	  // quadrature points of) the i-th shape function of current element.

	  // Once again, one can use "auto" instead of  "std::vector<FE_Function>".
	  // But in this case, it seems better to use a references, "auto&", because
	  // speed-up is improved (about 30%). Please test it.
	  auto& dx_phi = fe.get_dx_phi();
	  auto& dy_phi = fe.get_dy_phi();

	  Index ndofs = fe.get_ndofs();

	  // For each degree of freedom, i:
	  for (Index i = 0; i < ndofs; ++i)
	    {
	      // For each degree of freedom, j:
	      for (Index j = 0; j < ndofs; ++j)
		{
		  // Compute integral of product of x-derivatives
		  fe.integrate(dx_phi[i], dx_phi[j]);
		  // Compute integral of product of y-derivatives
		  fe.integrate(dy_phi[i], dy_phi[j]);
		}
	    }
	}
    }

  // Store final time and print elapsed time
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "Elapsed time: " <<
    std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() <<
    " miliseconds" << std::endl;
}
