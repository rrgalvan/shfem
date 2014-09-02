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

namespace shfem {

  /// Base class for Finite Element spaces on a given 2d mesh
  /// TODO: make mesh type a template argument
  class BaseFESpace {
  protected:
    BaseFESpace(Mesh const& m): mesh(&m) {}
    Mesh const* mesh;
  };

  /// Continous Galerkin ($P_k$--Lagrange) Finite Element spaces
  class CG_FESpace: 
    public BaseFESpace
  {
  public:
    CG_FESpace(Mesh const& m): BaseFESpace(m) {}
  };
  
}
#endif // FESPACE_HPP_
