// read_msh.cpp ---

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

#include <shfem.hpp>
#include <cstdlib>

int main()
{
  // Build an empty mesh object
  shfem::TriangleMesh mesh;

  // Try to read mesh contents from a .msh file (for information about the
  // structure of these files, see FreeFem++ documentation).
  try
    {
      mesh.read_file_msh("squared-mesh-2x2.msh");
    }
  catch (...)
    {
      std::cerr << "Error reading mesh file" << std::endl;
      exit(1);
    }

  // Show information about triangles and vertices
  std::cout << "Read mesh with " << mesh.get_ncel() << " triangles and "
	    << mesh.get_nver() << " vertices" << std::endl;

  // Print elements and coordinates of each vertex
  mesh.print();
}
