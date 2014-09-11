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

  /// Stores data for a concrete finite element
  class FiniteElementData;

  /// Base class for Finite Element spaces in a given 2d mesh
  /// TODO: make mesh type a template argument
  class BaseFESpace {
  protected:
    const Mesh* mesh;
    const BaseQuadRule* quadrature_rule;
  public:
    void set_mesh( const Mesh& m) { mesh = &m; }
    const Mesh& get_mesh() const { return *mesh; }
    void set_quadrature_rule(const BaseQuadRule& qr) { quadrature_rule = &qr; }
    const BaseQuadRule& get_quadrature_rule() const { return *quadrature_rule; }

    // const FiniteElementData& compute_data_for_element(index_t element_index) {
    // }
    
  };

  /// Continous Galerkin ($P_k$--Lagrange) Finite Element spaces
  class CG_FESpace: 
    public BaseFESpace
  {
  public:
    CG_FESpace(Mesh const& m) { set_mesh(m); }
  };
  
}

#endif // FESPACE_HPP_
