#ifndef QUADRATURE_HPP_
#define QUADRATURE_HPP_ 

#include "geometry_2d.hpp"

namespace shfem {
  

  class BaseQuadRule {
  };
  
  /// Available quadrature rules
  enum AvailableQuadRules { VerticesQR, MidpointsQR };

  /// Base class for Quadrature Rules on the reference element
  template <int QR>
  struct QuadRule : public BaseQuadRule {
    static std::vector<Point> nodes;
    static std::vector<real_t> weights;
  };

  // Vertices QR
  template<> std::vector<Point> QuadRule<VerticesQR>::nodes = 
   {Point(0.,0.), Point(0.,1.), Point(1.,0.)};
  template<> std::vector<real_t> QuadRule<VerticesQR>::weights = 
    {1./6, 1./6, 1./6};
}

#endif // QUADRATURE_HPP_ 
