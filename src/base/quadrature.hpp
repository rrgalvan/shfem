#ifndef QUADRATURE_HPP_
#define QUADRATURE_HPP_

#include "geometry_2d.hpp"

namespace shfem {

  class BaseQuadRule {
  };

  // A Quadfunction is defined by its values in the quadrature nodes
  typedef std::vector<Real> FE_Function;

  /// Available quadrature rules
  enum AvailableQuadRules { VerticesQR, MidpointsQR };

  namespace dim2 {

    /// Base class for Quadrature Rules on the reference element
    template <int QR>
    struct QuadRule : public BaseQuadRule {
      static std::vector<Point> nodes;
      static std::vector<Real> weights;

      /**
       * @brief Number of nodes and weights in current quadrature rule
       *
       * @return Number of nodes and weights
       */
      Index size() const { return weights.size();}

      /**
       * @brief Integrate a function (vector) on the reference elemenet
       * @param values vector representing values of the function on each
       *               quadrature point
       * @return Numerical approximation of the integral
       */
      Real integrate_on_ref_element(const std::vector<Real>& values) const;
    };

    // Nodes of the quadrature rule (on the reference element)
    template<> std::vector<Point> QuadRule<VerticesQR>::nodes =
      {Point(0.,0.), Point(1.,0.), Point(0.,1.)};
    // Weights of the quadrature rule
    template<> FE_Function QuadRule<VerticesQR>::weights =
      {1./6, 1./6, 1./6};

    template<int QR> Real
    QuadRule<QR>::integrate_on_ref_element(const FE_Function
					   & values) const {
      const int n=3;
      Real sum=0;
      for(int i=0; i<n;i++) {
	Real weight = this->weights[i];
	Real value = values[i];
	sum += weight * value;
      }
      return sum;
    }
  }
}
#endif // QUADRATURE_HPP_
