#ifndef QUADRATURE_HPP_
#define QUADRATURE_HPP_

#include "geometry_2d.hpp"

namespace shfem {

  /// Available quadrature rules
  enum AvailableQuadRules { VerticesQR, MidpointsQR };

  /// Base class for Quadrature Rules on the reference element
  template <int QR>
  struct QuadRule
  {
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

  // Nodes of the quadrature rule (on the reference element)
  template<> std::vector<Point> QuadRule<MidpointsQR>::nodes =
    {Point(0.5,0.), Point(0.5,0.5), Point(0.,0.5)};
  // Weights of the quadrature rule
  template<> FE_Function QuadRule<MidpointsQR>::weights =
    {1./6, 1./6, 1./6};

  template<int QR> Real
  QuadRule<QR>::integrate_on_ref_element(const FE_Function
					 & values) const {
    const int n=size();
    Real sum=0;
    for(int i=0; i<n;i++) {
      Real weight = weights[i];
      Real value = values[i];
      sum += weight * value;
    }
    return sum;
  }

  typedef QuadRule<VerticesQR> VerticesQuadRule;
}
#endif // QUADRATURE_HPP_
