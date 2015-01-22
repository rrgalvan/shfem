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

  /// Stores information for a concrete finite element
  class FiniteElement {
  private:
    const BaseQuadRule* _quadrature_rule;
  public:
    FiniteElement() :
      _quadrature_rule(NULL)
    {}
    /// Set quadrature rule
    /**
     * @brief Defines quadrature rule for current element
     * @param qr Quadrature rule
     */
    void set_quadrature_rule(const BaseQuadRule& qr) { _quadrature_rule = &qr; }

    /**
     * @brief Get quadrature rule attached to current element
     * @return Reference to quadrature rule (BaseQuadrule)
     */
    const BaseQuadRule& get_quadrature_rule() const { return *_quadrature_rule; }

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

  /// Base class for Finite Element spaces in a given 2d mesh
  //  Current implementation stores a vector of elements. It increases the need
  // for storage (and decreases computing requirements?)
  // - TODO: develop other implementations (not storing elements).
  // - TODO: make mesh type a template argument
  class BaseFESpace {
  public:
    typedef dim2::Mesh MeshType;
  protected:
    const MeshType* mesh;
    const BaseQuadRule* default_quadrature_rule;
    std::vector<FiniteElement> finite_elements; /// List of finite elements
  public:

    /**
     * Attach (a pointer to) a given mesh to current FE space and
     * adjusts consequently he size of finite_elements vector
     *
     * @param m Mesh to be attached
     */
    void set_mesh( const MeshType& m) {
      mesh = &m;
      finite_elements.resize(mesh->get_nelt());
    }

    /**
     * @brief Get Mesh
     * @return Mesh currently attached to this FESpace
     */
    const MeshType& get_mesh() const { return *mesh; }

    /// @brief Set quadrature rule that is used (by default) by all elements
    /// @param qr Quadrature rule
    void set_default_quadrature_rule(const BaseQuadRule& qr) {
      default_quadrature_rule = &qr; }

    /// @brief Get quadrature rule used (by default) by all elements
    /// @return Default quadrature rule
    const BaseQuadRule& get_default_quadrature_rule() const {
      return *default_quadrature_rule; }

    /**
     * @brief  Get a concrete finite element
     * @param element_index Index in current FESpace of the desired Element
     * @return Reference to the element
     */
    const FiniteElement& get_element(index_t element_index) const {
      return finite_elements[element_index];
    }
  };

  /// Continous Galerkin ($P_k$--Lagrange) Finite Element spaces
  class CG_FESpace:
    public BaseFESpace
  {
  public:
    typedef dim2::Mesh MeshType;
    CG_FESpace(MeshType const& m) { set_mesh(m); }
  };

}

#endif // FESPACE_HPP_
