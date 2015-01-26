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

  using FUNCTION_R2 = real_t(*)(real_t,real_t); // C++11 syntax

  struct ReferenceElement
  {
    typedef dim2::Point POINT;
    const std::vector<POINT> nodes;
    static real_t phi_0(real_t x, real_t y) {return 1-x-y;}
    static real_t phi_1(real_t x, real_t y) {return x;}
    static real_t phi_2(real_t x, real_t y) {return y;}
    const std::vector<FUNCTION_R2> basis_functions;
    ReferenceElement() :
      nodes({POINT(0.,0.), POINT(1.,0.), POINT(0.,1.)}),
      basis_functions({phi_0, phi_1, phi_2}) {}
  };

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
    typedef dim2::QuadRule<VerticesQR> QUADRULE;

    static ReferenceElement _reference_element;

  protected:
    const CELL* _cell;		/**< Topological information (global indices of vertex) */
    std::vector<POINT> _geometry; /**< Geometrical information (coord. of vertex) */
    const QUADRULE* _quadrature_rule; /**< Quadrature rule  */
    std::vector<FE_Function> _basis_functions; /**< Basis functions evaluated at quadrature points */

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
      size_t ndofs = get_ndofs();
      _geometry.resize(ndofs);
      _geometry[0] = mesh.get_vertex(idv0);
      _geometry[1] = mesh.get_vertex(idv1);
      _geometry[2] = mesh.get_vertex(idv2);

      // Compute element basis functions evaluated at the quadrature points.
      _basis_functions.resize(ndofs);
      for(size_t i=0; i<ndofs; ++i)
	{
	  // std::cout << "### idof=" << i << std::endl;
	  size_t nb_quad_nodes = _quadrature_rule->size();
	  _basis_functions[i].resize(nb_quad_nodes);
	  // // Get physical basis function (for i-th degree of freedom)
	  // const FE_Function& phi = _basis_functions[i];
	  // Get reference-element basis function (for i-th dof)
	  FUNCTION_R2 hat_phi = _reference_element.basis_functions[i];
	  // Loop on quadrature nodes
	  for(index_t inode=0; inode<nb_quad_nodes; ++inode)
	    {
	      // Get quadrature node (at reference triangle)
	      const POINT& hat_P = _quadrature_rule->nodes[inode];
	      // Apply reference basis function on that node
	      _basis_functions[i][inode] = hat_phi(hat_P.x, hat_P.y);
	    }
	  // for(index_t inode=0; inode<nb_quad_nodes; ++inode)
	  //   {
	  //     std::cout << "_basis_functions[" << i << "][" << inode <<  "] = "
	  // 		<< _basis_functions[i][inode] << std::endl;

	  //   }
	  // std::cout << std::endl;
	}
    }

    /**
     * @brief Get number of degrees of freedom for current element
     * @return Number of degrees of freedom
     */
    size_t get_ndofs() const {
      return NDOFS;
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
      real_t referenceX = point.x, referenceY = point.y;
      real_t x0 = get_vertex0().x, y0 = get_vertex0().y;
      real_t x1 = get_vertex1().x, y1 = get_vertex1().y;
      real_t x2 = get_vertex2().x, y2 = get_vertex2().y;
      float aux = 1-referenceX-referenceY;;
      return POINT(x0*aux + x1*referenceX + x2*referenceY,
		   y0*aux + y1*referenceX + y2*referenceY);
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
      real_t dx = x-x0, dy = y-y0;
      real_t referenceX = (dx*(y2-y0) + (x0-x2)*dy)*inv_det_B;
      real_t referenceY = (dx*(y0-y1) + (x1-x0)*dy)*inv_det_B;
      return POINT(referenceX, referenceY);
    }

    /**
     * @brief Get the element shape functions (evaluated at quadrature points)
     *
     * Each finite element contains a vector of basis functions (one
     * for each degree of freedom). Each basis function is defined as
     * a vector which stores the values of bais function at the
     * quadrature points.
     *
     * @return Vector of FEfunction (values in each quadrature  rule)
     */
    const std::vector<FE_Function>& get_basis_functions() const {
      return _basis_functions;
    }

    /**
     * @brief Return the basis function corresponding to a degree of freedom
     *
     * @param idof Index of the degree of freedom
     * @return FE_Function (vector of values at each quadrature point)
     */
    const FE_Function& get_basis_function(size_t idof) const {
      return _basis_functions[idof];
    }

    /**
     * @brief Get the x-dervative of element shape function
     * @return Vector of FEfunction (x-derivatives evaluated in each quadrature rule)
     */
    const std::vector<FE_Function>& get_dx_basis_functions() const;

    /**
     * @brief Get the y-dervative of element shape function
     * @return Vector of FEfunction (y-derivatives evaluated in each quadrature rule)
     */
    const std::vector<FE_Function>& get_dy_basis_functions() const;
  };

  ReferenceElement FiniteElement::_reference_element = ReferenceElement();
}

#endif // FINITE_ELEMENT_HPP_
