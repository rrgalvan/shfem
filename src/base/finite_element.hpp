// finite_element.hpp ---

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

#ifndef FINITE_ELEMENT_HPP_
#define FINITE_ELEMENT_HPP_

#include "geometry_2d.hpp"
#include "quadrature.hpp"

namespace shfem {

  /// @brief Stores information for a concrete finite element
  ///
  /// Information stored:
  ///   - Topology (global indices of each vertex)
  ///   - Geometry (coordinates of each vertex)
  ///   - Quadrature rule (to be used for numerical integration)
  ///
  /// TODO: define nb of dofs as a template parameter. Current
  /// implementation defines only 3 dof elements
  class FiniteElement {
  public:
    enum {NDOFS=3};
    typedef dim2::TriangleMesh MESH;
    typedef dim2::Point POINT;
    typedef MESH::CELL CELL;
    typedef std::vector<POINT> GEOMETRY;
    typedef BaseQuadRule QUADRULE;

  protected:
    const CELL* _cell;		/**< Topological information (global indices of vertex) */
    GEOMETRY _geometry;	/**< Geometrical information (coord. of vertex) */
    const QUADRULE* _quadrature_rule; /**< Quadrature rule  */

  public:
    /**
     * @brief Construct a void Finite element
     */
    FiniteElement() :
      _cell(NULL), _geometry(), _quadrature_rule(NULL) {}

    /**
     * @brief Delete data computed for current element ()
     */
    ~FiniteElement() {
      clear();
    }

    /**
     * @brief Compute element-specific data
     *
     * Computed data includes:
     *
     *   - Topological (global indices) and geometric (coordinates)
     *     information for vertices and also,
     *   - Evaluation of shape functions and of its derivatives on
     *     quadrature nodes)
     *
     * @param mesh Mesh for obtaining topological and geometric informtion
     * @param idx_cell Global index in the mesh of current cell (Element)
     */
    void reinit(const MESH& mesh, index_t idx_cell) {
      _cell = &mesh.get_cell(idx_cell);
      index_t idv0 = _cell->idv0; // Global index for vertex 1
      index_t idv1 = _cell->idv1; // Global index for vertex 2
      index_t idv2 = _cell->idv2; // Global index for vertex 3
      _geometry.resize(3);
      _geometry[0] = mesh.get_vertex(idv0);
      _geometry[1] = mesh.get_vertex(idv1);
      _geometry[2] = mesh.get_vertex(idv2);
    }

    /**
     * @brief Delete every data computed for current element (by function reinit())
     */
    void clear() {
      _cell = NULL;
      _geometry.clear();
      _quadrature_rule = NULL;
    }

    /**
     * @brief Defines underlying geometric cell for current element
     * @param cell Underlying cell to be attached
     */
    void set_cell(const CELL& cell) { _cell = &cell; }

    /**
     * @brief Returns the underlying topological cell for current element
     * @return Underlying cell
     */
    const CELL& get_cell() const { return *_cell; }

    /**
     * @brief Defines quadrature rule for current element
     * @param qr Quadrature rule
     */
    void set_quadrature_rule(const QUADRULE& qr) { _quadrature_rule = &qr; }

    /**
     * @brief Get quadrature rule attached to current element
     * @return Reference to quadrature rule (BaseQuadrule)
     */
    const QUADRULE& get_quadrature_rule() const { return *_quadrature_rule; }

    const POINT& get_vertex0() const { return _geometry[0]; }
    const POINT& get_vertex1() const { return _geometry[1]; }
    const POINT& get_vertex2() const { return _geometry[2]; }
    const POINT& get_vertex(index_t i) const { return _geometry[i]; }


    // Returns $F_T(\hat P)$, where $T$ is this triangle, $\hat P \in T$ and
    // $F_T$ is the affine transformation from the reference cell to $T$.
    // Point affine_transform(const Triangle& T, const&  hatP) const; // UNIMPLEMENTED

    /**
     * @brief Determinant of Jacobian of map $F_K$
     *
     * Calculate the determinant of the Jacobian of the
     * affine transformation $F_K$ for a given cell K
     *
     * @return Determinant of Jacobian
     */
    real_t det_jacobian_affine_map() const {
      real_t x0 = get_vertex0().x, y0 = get_vertex0().y;
      real_t x1 = get_vertex1().x, y1 = get_vertex1().y;
      real_t x2 = get_vertex2().x, y2 = get_vertex2().y;
      return (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0);
    }

    /**
     * Calculate the area of a (triangle) cell
     *
     * @return Area of cell
     */
    real_t area() const {
      static const real_t area_of_reference_cell = 0.5;
      return area_of_reference_cell * std::abs(det_jacobian_affine_map());
    }

    /**
     * @brief Apply the affine map to a point located in the reference cell
     *
     * If $K$ is the cell attached to current element, the affine map
     * is defined as the application $T_K: \overline K \to K$ (where
     * $\overline K$ is a reference polyhedron) such that
     * $$T_K(\overline x) = B_K \overline x + b_K, $$ $B_K$ being a
     * nonsingular matrix.
     *
     * @param point Point located in the reference cell
     * @return Point resulting from application of the affine map
     */
    POINT apply_affine_map(const POINT& point) {
      real_t hatx = point.x, haty = point.y;
      real_t x0 = get_vertex0().x, y0 = get_vertex0().y;
      real_t x1 = get_vertex1().x, y1 = get_vertex1().y;
      real_t x2 = get_vertex2().x, y2 = get_vertex2().y;
      float aux = 1-hatx-haty;;
      return POINT(x0*aux + x1*hatx + x2*haty,
		   y0*aux + y1*hatx + y2*haty);
    }

    /**
     * @brief Apply inverse affine map to a point located in the physical cell
     *
     * If $K$ is the cell attached to current element, the affine map
     * is defined as the application $T^{-1}_K: K \to \overline K$
     * (where $\overline K$ is a reference polyhedron) such that
     * $$T_K(\overline x) = B_K \overline x + b_K, $$ $B_K$ being a
     * nonsingular matrix.
     *
     * @param point Point located in the physical cell
     * @return Point resulting from application of the inverse affine map
     */
    POINT apply_inv_affine_map(const POINT& point)
    {
      real_t x = point.x, y = point.y;
      real_t x0 = get_vertex0().x, y0 = get_vertex0().y;
      real_t x1 = get_vertex1().x, y1 = get_vertex1().y;
      real_t x2 = get_vertex2().x, y2 = get_vertex2().y;
      real_t inv_det_B = 1.0/((x1-x0)*(y2-y0) - (y1-y0)*(x2-x0));
      real_t dx = x-x0, dy=y-y0;
      real_t hatx = (dx*(y2-y0) + (x0-x2)*dy)*inv_det_B;
      real_t haty = (dx*(y0-y1) + (x1-x0)*dy)*inv_det_B;
      return POINT(hatx, haty);
    }

    /**
     * @brief Get the element shape functions (evaluated at quadrature points)
     * @return Vector of FEfunction (values in each quadrature rule)
     */
    const std::vector<QuadFunction>& get_shape_functions() const;

    /**
     * @brief Get the x-dervative of element shape function
     * @return Vector of FEfunction (x-derivatives evaluated in each quadrature rule)
     */
    const std::vector<QuadFunction>& get_dx_shape_functions() const;

    /**
     * @brief Get the y-dervative of element shape function
     * @return Vector of FEfunction (y-derivatives evaluated in each quadrature rule)
     */
    const std::vector<QuadFunction>& get_dy_shape_functions() const;
  };

}

#endif // FINITE_ELEMENT_HPP_
