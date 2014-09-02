// geometry_2d.hpp --- 

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

#ifndef GEOMETRY_2D_HPP_
#define GEOMETRY_2D_HPP_

#include "hpfem_base.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

namespace hpfem {
  
  //,---------
  //| Point 2D
  //`---------
  struct Point {
    Point() {}
    Point(real_t xx, real_t yy) : x(xx), y(yy) {}
    real_t x;
    real_t y;
    void print() const;
  };

  void Point::print() const {
    std::cout << "(" << this->x << "," << this->y << ")" << std::endl;
  }

  //,---------
  //| Triangle
  //`---------
  struct Triangle {
    Triangle() {}
    Triangle(index_t i1, index_t i2, index_t i3): idv1(i1), idv2(i2), idv3(i3) {}
    // Global index of vertex 1
    index_t idv1; 
    // Global index of vertex 2
    index_t idv2;
    // Global index of vertex 3
    index_t idv3;
  };


  //,-----
  //| Mesh
  //`-----
  class Mesh {
    std::vector<Point> vertices;
    std::vector<Triangle> triangles;
    

  public:
    /// Return number of vertices
    index_t get_nv() const { return vertices.size(); }
    /// Return number of triangles
    index_t get_nt() const { return triangles.size(); }
    /// Read mesh from a medit .msh file (e.g. wrote by FreeFem++)
    void read_file_msh(const char* filename);
    /// Print mesh contents
    void print() const;    

    /// Returns $F_T(\hat P)$, where $T$ is this triangle, $\hat P \in T$ and
    /// $F_T$ is the affine transformation from the reference element to $T$.
    // Point affine_transform(const Triangle& T, const&  hatP) const; // UNIMPLEMENTED

    /// Returns the determinant of the Jacobian of the affine transformation $F_T$
    real_t det_J_affine_transform(index_t triangle_id) const {
      const Triangle& T = triangles[triangle_id];
      real_t x1 = vertices[T.idv1].x;
      real_t y1 = vertices[T.idv1].y;
      real_t x2 = vertices[T.idv2].x;
      real_t y2 = vertices[T.idv2].y;
      real_t x3 = vertices[T.idv3].x;
      real_t y3 = vertices[T.idv3].y;
      return (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1);
    }
    
    real_t area(index_t triangle_id) const {
      return 2*std::abs(det_J_affine_transform(triangle_id));
    }
  };
  
  void Mesh::print() const {
    for(index_t i=0; i<get_nt(); ++i) {
      std::cout << "Triangle " << i << ":" << std::endl;
      vertices[ triangles[i].idv1 ].print();
      vertices[ triangles[i].idv2 ].print();
      vertices[ triangles[i].idv3 ].print();
    }
  }

  void Mesh::read_file_msh(const char* filename) {
    std::fstream meshfile(filename);
  
    int nv, nt, ne;
    meshfile >> nv >> nt >> ne;
    this->vertices = std::vector<Point>(nv);
    this->triangles = std::vector<Triangle>(nt);

    for(int i=0; i<nv; i++) {
      real_t x, y, label;
      meshfile >> x >> y >> label;
      Point p(x, y);
      vertices[i]=p;  
    }
  
    for(int i=0; i<nt; i++) {
      unsigned int id1, id2, id3, zero;
      meshfile >> id1 >> id2 >> id3 >> zero;
      triangles[i] = Triangle(id1-1, id2-1, id3-1); // 1 index -> 0 index
    }
  }

}

#endif // GEOMETRY_2D_HPP_
