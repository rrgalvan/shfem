// fespace.hpp ---

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

#ifndef FESPACE_HPP_
#define FESPACE_HPP_

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
    typedef MESH::CELL CELL;
    typedef dim2::TriangleGeometry GEOMETRY;
    typedef BaseQuadRule QUADRULE;

  protected:
    const CELL* _cell;		/**< Topological information (global indices of vertex) */
    const GEOMETRY* _geometry;	/**< Geometrical information (coord. of vertex) */
    const QUADRULE* _quadrature_rule; /**< Quadrature rule  */

  public:
    /**
     * @brief Construct a void Finite element
     */
    FiniteElement() :
      _cell(NULL), _geometry(NULL), _quadrature_rule(NULL) {}

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
    void reinit(const MESH& mesh, index_t idx_cell){
      _cell = &mesh.get_cell(idx_cell);
      index_t idv1 = _cell->idv1; // Global index for vertex 1
      index_t idv2 = _cell->idv2; // Global index for vertex 2
      index_t idv3 = _cell->idv3; // Global index for vertex 3
      const dim2::Point& v1 = mesh.get_vertex(idv1);
      const dim2::Point& v2 = mesh.get_vertex(idv2);
      const dim2::Point& v3 = mesh.get_vertex(idv3);
      _geometry = new GEOMETRY(v1, v2, v3);
    }

    /**
     * @brief Delete every data computed for current element (by function reinit())
     */
    void clear() {
      _cell = NULL;
      if (_geometry != NULL) { delete _geometry; _geometry = NULL; }
      _quadrature_rule = NULL;
    }

    /**
     * @brief Defines underlying geometric cell for current element
     * @param cell Underlying cell to be attached
     */
    void set_cell(const CELL& cell) { _cell = &cell; }

    /**
     * @brief Returns the underlying geometric cell for current element
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


  // /**
  //  * @brief Base class for Finite Element spaces in a given 2d mesh
  //  *
  //  * Current implementation stores a vector of elements. It increases
  //  * the need for storage (and decreases computing requirements?)
  //  *
  //  * - TODO: develop other implementations (not storing elements).
  //  * - TODO: make mesh type a template argument
  //  */
  // class BaseFESpace {
  //   // public:
  //   //   typedef dim2::TriangleMesh MeshType;
  //   //   typedef FiniteElement ELEMENT;
  // protected:
  //   const QUADRULE* default_quadrature_rule;
  //   std::vector<ELEMENT> finite_elements; /// List of finite elements
  // public:
  //   BaseFESpace(): mesh(NULL), default_quadrature_rule(NULL) {}

  //   /**
  //    * Attach (a pointer to) a given mesh to current FE space.

  //    * Also adjust consequently he size of finite_elements vector (in
  //    * default implementation, it is defined as the number of cells in
  //    * the atached mesh)
  //    *
  //    * @param m Mesh to be attached
  //    */
  //   void set_mesh( const MESH& m) {
  //     mesh = &m;
  //     int n = mesh->get_ncel();
  //     finite_elements.resize(n);
  //     for(int i=0; i<n; i++) {
  // 	ELEMENT* fe = &finite_elements[i];
  // 	fe->set_cell(mesh->get_cell(i));
  // 	fe->set_quadrature_rule(*default_quadrature_rule);
  //     }
  //   }

  //   /**
  //    * @brief Get Mesh
  //    * @return Mesh currently attached to this FESpace
  //    */
  //   const MeshType& get_mesh() const { return *mesh; }

  //   /// @brief Set quadrature rule that is used (by default) by all elements
  //   /// @param qr Quadrature rule
  //   void set_default_quadrature_rule(const QUADRULE& qr) {
  //     default_quadrature_rule = &qr; }

  //   /// @brief Get quadrature rule used (by default) by all elements
  //   /// @return Default quadrature rule
  //   const QUADRULE& get_default_quadrature_rule() const {
  //     return *default_quadrature_rule; }

  //   /**
  //    * @brief Get the number of elements in current FE space.
  //    *
  //    * For the default implementation, the number of elements is
  //    * defined as the size of the finite_elements list
  //    *
  //    * @return Number of elements in current FE space
  //    */
  //   size_t get_nelt() const {
  //     return this->finite_elements.size();
  //   }

  //   /**
  //    * @brief Get a concrete finite element
  //    * @param element_index Index in current FESpace of the desired Element
  //    * @return Reference to the element
  //    */
  //   const ELEMENT& get_element(index_t element_index) const {
  //     return finite_elements[element_index];
  //   }
  // };

  // /// Continous Galerkin ($P_k$--Lagrange) Finite Element spaces
  // class CG_FESpace:
  //   public BaseFESpace
  // {
  // public:
  //   typedef dim2::TriangleMesh MeshType;
  //   CG_FESpace() {}
  //   // CG_FESpace(MeshType const& m) { set_mesh(m); }
  // };

}

#endif // FESPACE_HPP_
