// finite_element.hpp ---

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

#ifndef FINITE_ELEMENT_HPP_
#define FINITE_ELEMENT_HPP_

#include "geometry_2d.hpp"
#include "quadrature.hpp"

#include <cassert> // For assert()

namespace shfem {

  using FUNCTION_R2R1 = Real(*)(Real,Real); // Function: R^2->R (C++11 syntax)

  struct ReferenceElement
  {
    typedef Point POINT;
    const std::vector<POINT> nodes;
    // Basis functions
    static Real phi_0(Real x, Real y) {return 1-x-y;}
    static Real phi_1(Real x, Real y) {return x;}
    static Real phi_2(Real x, Real y) {return y;}
    // Derivatives of basis functions
    static Real dx_phi_0(Real x, Real y) {return -1;}
    static Real dx_phi_1(Real x, Real y) {return x;}
    static Real dx_phi_2(Real x, Real y) {return 0;}
    static Real dy_phi_0(Real x, Real y) {return -1;}
    static Real dy_phi_1(Real x, Real y) {return 0;}
    static Real dy_phi_2(Real x, Real y) {return 1;}

    const std::vector<FUNCTION_R2R1> basis_functions;
    const std::vector<FUNCTION_R2R1> dx_phi;
    const std::vector<FUNCTION_R2R1> dy_phi;

    ReferenceElement() :
      nodes({POINT(0.,0.), POINT(1.,0.), POINT(0.,1.)}),
      basis_functions({phi_0, phi_1, phi_2}),
      dx_phi({dx_phi_0, dx_phi_1, dx_phi_2}),
      dy_phi({dy_phi_0, dy_phi_1, dy_phi_2}) {}
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
    typedef TriangleMesh MESH;
    typedef Point POINT;
    typedef MESH::CELL CELL;
    typedef QuadRule<VerticesQR> QUADRULE;

    static ReferenceElement _reference_element;

  protected:
    const CELL* _cell;	/**< Topological information (global indices of vertex) */
    std::vector<POINT> _geometry; /**< Geometrical information (coord. of vertex) */
    const QUADRULE* _quadrature_rule; /**< Quadrature rule  */
    std::vector<FE_Function> _phi; /**< Basis functions evaluated on points */
    std::vector<FE_Function> _dx_phi; /**< x-derivative of basis functions on quadrature points */
    std::vector<FE_Function> _dy_phi; /**< y_derivative of basis functions on quadrature points */

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
     *   - Evaluation of shape functions (and of its derivatives) on
     *     quadrature nodes
     *
     * @param mesh Mesh for obtaining topological and geometric informtion
     * @param idx_cell Global index in the mesh of current cell (Element)
     */
    void reinit(const MESH& mesh, Index idx_cell) {
      //,--------------------------------------------------
      //| Compute global index and coordinates for vertices
      //`--------------------------------------------------
      _cell = &mesh.get_cell(idx_cell);
      Index idv0 = _cell->idv0; // Global index for vertex 0
      Index idv1 = _cell->idv1; // Global index for vertex 1
      Index idv2 = _cell->idv2; // Global index for vertex 2

      Index ndofs = get_ndofs();
      _geometry.resize(ndofs);
      _geometry[0] = mesh.get_vertex(idv0); // Coord. of vertex 0
      _geometry[1] = mesh.get_vertex(idv1); // Coord. of vertex 2
      _geometry[2] = mesh.get_vertex(idv2); // Coord. of vertex 2

      //,------------------------------------------------------------------
      //| Evaluate basis functions (and its derivatives) on quadrature nodes
      //`------------------------------------------------------------------
      _phi.resize(ndofs);
      _dx_phi.resize(ndofs);
      _dy_phi.resize(ndofs);
      for(Index i=0; i<ndofs; ++i)
	{
	  // std::cout << "### idof=" << i << std::endl;
	  Index nb_quad_nodes = _quadrature_rule->size();
	  _phi[i].resize(nb_quad_nodes);
	  _dx_phi[i].resize(nb_quad_nodes);
	  _dy_phi[i].resize(nb_quad_nodes);
	  // Define i-th basis function (on the reference-element)
	  FUNCTION_R2R1 hat_phi_i = _reference_element.basis_functions[i];
	  FUNCTION_R2R1 dx_hat_phi_i = _reference_element.dx_phi[i];
	  FUNCTION_R2R1 dy_hat_phi_i = _reference_element.dy_phi[i];
	  // Loop on quadrature nodes
	  for(Index j=0; j<nb_quad_nodes; ++j)
	    {
	      // Get quadrature node (at reference triangle)
	      const POINT& hat_P_j = _quadrature_rule->nodes[j];
	      // Apply reference basis function on that node
	      _phi[i][j] = hat_phi_i(hat_P_j.x, hat_P_j.y);
	      _dx_phi[i][j] = dx_hat_phi_i(hat_P_j.x, hat_P_j.y);
	      _dy_phi[i][j] = dy_hat_phi_i(hat_P_j.x, hat_P_j.y);
	    }
	}
    }

    /**
     * @brief Get number of degrees of freedom for current element
     * @return Number of degrees of freedom
     */
    Index get_ndofs() const {
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
    const POINT& get_vertex(Index i) const { return _geometry[i]; }

    /**
     * @brief Get the element shape functions (evaluated on quadrature points)
     *
     * Each finite element contains a vector of basis functions (one
     * for each degree of freedom). Each basis function is defined as
     * a FE_Function which stores the values of basis function at the
     * quadrature points.
     *
     * @return Vector of FE_function (values on each quadrature  rule)
     */
    const std::vector<FE_Function>& get_phi() const {
      return _phi;
    }

    /**
     * @brief Get x-derivative of element shape functions (on quadrature points)
     *
     * Each derivative function is defined as a FE_Function, which stores
     * the values of the derivative at the quadrature points.
     *
     * @return Vector of FE_function (values of derivative on each quadrature rule)
     */
    const std::vector<FE_Function>& get_dx_phi() const {
      return _dx_phi;
    }

    /**
     * @brief Get y-derivative of element shape functions (on quadrature points)
     *
     * Each derivative function is defined as a FE_Function, which stores
     * the values of the derivative at the quadrature points.
     *
     * @return Vector of FE_function (values of derivative on each quadrature rule)
     */
    const std::vector<FE_Function>& get_dy_phi() const {
      return _dy_phi;
    }

    /**
     * @brief Return the basis function corresponding to a degree of freedom
     *
     * @param idof Index of the degree of freedom
     * @return FE_Function (vector of values at each quadrature point)
     */
    const FE_Function& get_phi(Index idof) const {
      return _phi[idof];
    }

    /**
     * @brief Get the x-dervative of the shape function corresponding to a dof
     * @return Vector of FEfunction (x-derivatives evaluated in each quadrature rule)
     */
    const FE_Function& get_dx_phi(Index idof) const {
      return _dx_phi[idof];
    }

    /**
     * @brief Get the y-dervative of the shape function corresponding to a dof
     * @return Vector of FEfunction (y-derivatives evaluated in each quadrature rule)
     */
    const FE_Function& get_dy_phi(Index idof) const {
      return _dy_phi[idof];
    }

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
    Real det_jacobian_affine_map() const {
      Real x0 = get_vertex0().x, y0 = get_vertex0().y;
      Real x1 = get_vertex1().x, y1 = get_vertex1().y;
      Real x2 = get_vertex2().x, y2 = get_vertex2().y;
      return (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0);
    }

    /**
     * Calculate the area of a (triangle) cell
     *
     * @return Area of cell
     */
    Real area() const {
      static const Real area_of_reference_cell = 0.5;
      return area_of_reference_cell * std::abs(det_jacobian_affine_map());
    }

    /**
     * @brief Computes the integral of f1*f2 in current element
     *
     * Integral is aproximated as
     *  $$\int_{K} f_1\dot f_2 = \sum{i=0}^{N} w_i f_1(p_i) f_2(p_i),$$
     * where $K$ is current element, $N+1$ is the number of quadrature points, and
     * $w_i$, $p_i$ are wheights and points, respectively.
     *
     * @param f1 FE_Function (vector of reals) evaluated on quadrature points
     * @param f2 FE_Function (vector of reals) evaluated on quadrature points
     *
     * @return Numerical approximation of $\int_{K} f_1\dot f_2$
     */
    Real integrate(const FE_Function& f1, const FE_Function& f2) const
    {
      Real sum = 0.;
      Real abs_det_J = std::abs(det_jacobian_affine_map());
      Index nb_points = _quadrature_rule->size();
      const std::vector<Real>& w = _quadrature_rule->weights;
      for (Index i = 0; i < nb_points; ++i)
	{
	  sum += w[i]*f1[i]*f2[i];
	}
      return abs_det_J*sum;
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
    POINT apply_affine_map(const POINT& point) const {
      Real referenceX = point.x, referenceY = point.y;
      Real x0 = get_vertex0().x, y0 = get_vertex0().y;
      Real x1 = get_vertex1().x, y1 = get_vertex1().y;
      Real x2 = get_vertex2().x, y2 = get_vertex2().y;
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
    POINT apply_inv_affine_map(const POINT& point) const
    {
      Real x = point.x, y = point.y;
      Real x0 = get_vertex0().x, y0 = get_vertex0().y;
      Real x1 = get_vertex1().x, y1 = get_vertex1().y;
      Real x2 = get_vertex2().x, y2 = get_vertex2().y;
      Real inv_det_B = 1.0/((x1-x0)*(y2-y0) - (y1-y0)*(x2-x0));
      Real dx = x-x0, dy = y-y0;
      Real referenceX = (dx*(y2-y0) + (x0-x2)*dy)*inv_det_B;
      Real referenceY = (dx*(y0-y1) + (x1-x0)*dy)*inv_det_B;
      return POINT(referenceX, referenceY);
    }

  };


  /**
   * @brief Espace of P1 continous finite elements
   */
  class P1_FE_Space {
  public:
    typedef typename FiniteElement::MESH MESH;
    typedef typename FiniteElement::QUADRULE QUADRULE;

    P1_FE_Space(MESH mesh, QUADRULE quad_rule):
      _mesh(mesh), _default_quadrature_rule(quad_rule) {}

    const FiniteElement get_element(Index r) const {
      FiniteElement fe;
      fe.reinit(_mesh, r);
      // WARNING! USE MOVE SEMANTICS.
      return fe;
    }

    /**
     * @brief Get number of degrees of freedom in current space
     * @return Number of degrees of freedom
     */
    Index get_ndofs() const {
      return _mesh.get_nver(); // P1 dofs match mesh vertices
    }

    // Assemble local matrix "Mlocal" in global matrix "M"
    template <class Matrix> void assemble_matrix(const Matrix& Mlocal, Matrix& M) const
    {
      // For each cell, r:
      for (Index r=0; r<_mesh.get_ncel(); ++r) {
	//,--------------------------------------------------
	//| Compute global index and coordinates for vertices
	//`--------------------------------------------------
	auto _cell = _mesh.get_cell(r);
	Index idv0 = _cell.idv0; // Global index for vertex 0
	Index idv1 = _cell.idv1; // Global index for vertex 1
	Index idv2 = _cell.idv2; // Global index for vertex 2
	const Index index_map[3] = {idv0, idv1, idv2};
	const Index ndofs = _cell.get_nver();
	assert(ndofs==3);
	for (Index i = 0; i < ndofs; ++i)
	  {
	    for (Index j = 0; j < ndofs; ++j)
	      {
		// Here we are assuming that class Matrix implements
		// operator()(int i,int j) for accessing to element (i,j)
		M(index_map[i], index_map[j]) += Mlocal(i,j);
	      }
	  }
      }
    }

  private:
    // We do not need storing pointers or refferences because we
    // assume that MESH and QUADRULE uses C++11 move semantics.
    MESH _mesh;
    QUADRULE _default_quadrature_rule;
  };

  //
  // Define static member FiniteElement::_reference_element
  //
  ReferenceElement FiniteElement::_reference_element = ReferenceElement();
}

#endif // FINITE_ELEMENT_HPP_
