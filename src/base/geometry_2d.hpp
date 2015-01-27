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

#include "shfem_base.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

namespace shfem {

  class BaseMesh {};

  namespace dim2 {

    //,---------
    //| Point 2D
    //`---------
    struct Point {
      Point() {}
      Point(Real xx, Real yy) : x(xx), y(yy) {}
      Real x;
      Real y;
      void print() const;
    };

    /**
     * Print coordinates of current point
     */
    void Point::print() const {
      std::cout << "(" << this->x << "," << this->y << ")" << std::endl;
    }


    // /**
    //  * @brief Topological information for a 2d triangle.
    //  *
    //  * Class wich stores three indices corresponding to the vertices
    //  * of a 2d tirnagle
    //  */
    // struct Triangle {
    //   Triangle() {}
    //   Triangle(Index i1, Index i2, Index i3): idv1(i1), idv2(i2), idv3(i3) {}
    //   // Global index of vertex 1
    //   Index idv1;
    //   // Global index of vertex 2
    //   Index idv2;
    //   // Global index of vertex 3
    //   Index idv3;
    // };

    /**
     * @brief Topological information for a 2d triangle.
     *
     * Class wich stores three indices corresponding to the global
     * information about vertices of a 2d trangle
     */
    struct Triangle {
      Triangle() {}
      Triangle(Index i0, Index i1, Index i2): idv0(i0), idv1(i1), idv2(i2) {}
      // Global index of vertex 0
      Index idv0;
      // Global index of vertex 1
      Index idv1;
      // Global index of vertex 2
      Index idv2;
    };

    /**
     * @brief Geometric information for a topological object
     *
     * Class wich stores referencies to three points, storing the
     * coordinates of vertices of topological object (for instance, a
     * triangle)
     */
    struct TriangleGeometry {
      const Point& vertex0;
      const Point& vertex1;
      const Point& vertex2;

      /**
       * @brief Buld geometry, storing references to three Points
       * @param p0 Coordinates of first point
       * @param p1 Coordinates of second point
       * @param p2 Coordinates of third point
       */
      TriangleGeometry(const Point& p0, const Point& p1, const Point& p2):
	vertex0(p0), vertex1(p1), vertex2(p2) {}
    };

    /**
     * @brief Mesh composed by 2d triangles
     */
    class TriangleMesh: public BaseMesh {
    public:
      typedef Triangle CELL;
    private:
      std::vector<Point> vertices;
      std::vector<CELL> cells;
    public:
      /// @brief Read number of vertices in current mesh
      /// @return number of vertices
      Index get_nver() const { return vertices.size(); }

      /// @brief Read number of cells in current mesh
      /// @return number of cells
      Index get_ncel() const { return cells.size(); }

      /// @brief Get a concrete cell contained in current mesh
      /// @return (Reference to) a cell
      const CELL& get_cell(Index i) const { return cells[i];}

      /// @brief Get a concrete vertex contained in current mesh
      /// @return (Reference to) a Point
      const Point& get_vertex(Index i) const { return vertices[i];}

      /// @brief Read mesh from a medit .msh file (e.g. wrote by FreeFem++)
      /// @param filename Name of a file containing the mesh
      void read_file_msh(const char* filename);

      /// @brief Print mesh contents
      void print() const;

    };

    /**
     * Print information about current mesh
     *
     * Displayed information: identifier for each cell
     * and vertices of each cell
     */
    void TriangleMesh::print() const {
      for(Index i=0; i<get_ncel(); ++i) {
	const Triangle& T = cells[i];
	std::cout << "Cell " << i << ":" << std::endl;
	vertices[ T.idv0 ].print();
	vertices[ T.idv1 ].print();
	vertices[ T.idv2 ].print();
      }
    }

    /**
     * Read mesh (msh format) from a file
     *
     * Try to read mesh contents from a .msh file (for information about the
     * structure of these files, see FreeFem++ documentation).
     *
     * @param filename
     */
    void TriangleMesh::read_file_msh(const char* filename) {
      std::fstream meshfile(filename);

      int nver, ncel, nedg;
      meshfile >> nver >> ncel >> nedg;
      this->vertices = std::vector<Point>(nver);
      this->cells = std::vector<CELL>(ncel);

      for(int i=0; i<nver; i++) {
	Real x, y, label;
	meshfile >> x >> y >> label;
	Point p(x, y);
	vertices[i]=p;
      }

      for(int i=0; i<ncel; i++) {
	unsigned int id0, id1, id2, zero;
	meshfile >> id0 >> id1 >> id2 >> zero;
	cells[i] = CELL(id0-1, id1-1, id2-1); // 1 index -> 0 index
      }
    }

  }

  typedef dim2::Point Point2D;
}

#endif // GEOMETRY_2D_HPP_
