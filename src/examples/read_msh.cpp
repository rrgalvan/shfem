// read_msh.cpp --- 

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

#include <hpfem.hpp>

int main()
{
  // Build an empty mesh object
  hpfem::Mesh m;
  
  // Try to read mesh contents from a .msh file (for information about the
  // structure of these files, see FreeFem++ documentation).
  try{
    m.read_file_msh("squared-mesh-2x2.msh");
  } catch (...) {
    std::cerr << "Error reading mesh file" << std::endl;
  }

  // Show information about triangles and vertices
  std::cout << "Read mesh with " << m.get_nt() << " triangles and "
	    << m.get_nv() << " vertices" << std::endl;

  // Print elements and coordinates of each vertex
  m.print();
  
  // Calculate area of each triangle 
  for(unsigned i=0; i<m.get_nt(); ++i) {
    std::cout << "Area of triangle " << i << ": " << m.area(i) << std::endl;
  }
}
