// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <Eigen/Core>
#include <dune/common/fvector.hh>
#include <dune/localfefunctions/impl/standardLocalFunction.hh>
#include <dune/common/overloadset.hh>
#include <ranges>
namespace Ikarus {



struct DefaultMembraneStrain {

  template<typename Geometry>
  void pre( const Geometry &geo,
            const auto &uFunction){}

  template< typename Geometry>
  auto value(const Dune::FieldVector<double, 2> &gpPos,
             const Geometry &geo,
             const auto &uFunction) const -> Eigen::Vector3<typename std::remove_cvref_t<decltype(uFunction)>::ctype>{
    using ScalarType = typename std::remove_cvref_t<decltype(uFunction)>::ctype;
    Eigen::Vector3<ScalarType> epsV;
    const auto J = Dune::toEigen(geo.jacobianTransposed(gpPos));
    using namespace Dune;
    using namespace Dune::DerivativeDirections;
    const Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(
        uFunction.evaluateDerivative(gpPos, //Here the gpIndex could be passed
                                     Dune::wrt(spatialAll),
                                     Dune::on(Dune::DerivativeDirections::referenceElement)));
    const Eigen::Matrix<ScalarType, 2, 3> j = J + gradu.transpose();

    epsV << J.row(0).dot(gradu.col(0)) + 0.5*gradu.col(0).squaredNorm(),
        J.row(1).dot(gradu.col(1)) + 0.5*gradu.col(1).squaredNorm(), j.row(0).dot(j.row(1));
    return epsV;
  }
  template<typename Geometry,typename ScalarType>
  auto derivative(const Dune::FieldVector<double, 2> &gpPos,const Eigen::Matrix<ScalarType, 2, 3> &jcur, const auto &dNAtGp, const Geometry& geo,const auto& uFunction, const auto& localBasis,
                  const int node)  const{
    Eigen::Matrix<ScalarType, 3, 3> bop;
    bop.row(0) = jcur.row(0)*dNAtGp(node, 0);
    bop.row(1) = jcur.row(1)*dNAtGp(node, 1);
    bop.row(2) = jcur.row(0)*dNAtGp(node, 1) + jcur.row(1)*dNAtGp(node, 0);

    return bop;
  }
  template<typename Geometry,typename ScalarType>
  auto secondDerivative(const Dune::FieldVector<double, 2> &gpPos,const auto &dNAtGp, const Geometry& geo,const auto& uFunction, const auto& localBasis,
                        const Eigen::Vector3<ScalarType> &S, int I, int J)const {
    const auto &dN1i = dNAtGp(I, 0);
    const auto &dN1j = dNAtGp(J, 0);
    const auto &dN2i = dNAtGp(I, 1);
    const auto &dN2j = dNAtGp(J, 1);
    const ScalarType NS = dN1i*dN1j*S[0] + dN2i*dN2j*S[1] + (dN1i*dN2j + dN2i*dN1j)*S[2];
    Eigen::Matrix<ScalarType, 3, 3> kg = Eigen::Matrix<double, 3, 3>::Identity()*NS;
    return kg;
  }

};


struct CASAnsatzFunction
{
  Dune::LagrangeCubeLocalFiniteElement<double, double, 2, 1> q1lfem2D;
  mutable std::vector<double> out;

  auto positions(std::vector<Dune::FieldVector<double, 2>>& lagrangePoints)
  {
    lagrangePoints.resize(q1lfem2D.size());
    for (int i = 0; i < 2; i++) {
      auto ithCoord = [&i](const Dune::FieldVector<double, 2>& x) { return x[i]; };
      q1lfem2D.localInterpolation().interpolate(ithCoord, out);
      for (std::size_t jI = 0; jI < out.size(); jI++)
        lagrangePoints[jI][i] = out[jI];
    }
  }

  void evaluateFunction(const Dune::FieldVector<double,2>& gpPos,std::vector<Dune::FieldVector<double, 1>>& N) const
  {
    q1lfem2D.localBasis().evaluateFunction(gpPos, N);
  }
};

struct CASAnsatzFunctionANS
{
  Dune::LagrangeCubeLocalFiniteElement<double, double, 2, 1> q1lfem2D;
  mutable std::vector<double> out;

  auto positions(std::vector<Dune::FieldVector<double, 2>>& lagrangePoints)
  {
    lagrangePoints.resize(4);
    lagrangePoints[0] = {0, 0.5};
    lagrangePoints[1] = {0.5, 0};
    lagrangePoints[2] = {1, 0.5};
    lagrangePoints[3] = {0.5, 1};
  }

  void evaluateFunction(const Dune::FieldVector<double,2>& gpPos,std::vector<Dune::FieldVector<double, 1>>& NANS) const
  {
    NANS.resize(4);
    NANS[0] =  (1 - gpPos[1]);
    NANS[1] = (1 - gpPos[0]);
    NANS[2] = gpPos[1];
    NANS[3] = gpPos[0];
  }

};

template<typename Ansatzfunction>
struct CASMembraneStrain
{
  DefaultMembraneStrain defaultMembraneStrain;
  mutable std::vector<Dune::FieldVector<double, 2>> lagrangePoints;
  Ansatzfunction cASAnsatzFunction;


  template<typename T>
  using Vec = std::vector<Eigen::Vector<T, 3>>;
  std::variant<Vec<double>,Vec<autodiff::dual>,Vec<autodiff::dual2nd>> membraneStrainsAtVertices;
  mutable std::vector<Dune::FieldVector<double, 1>> NANS;

  CASMembraneStrain()
  {
    cASAnsatzFunction.positions(lagrangePoints);
}
  template<typename Geometry>
  void pre( const Geometry &geo,
           const auto &uFunction) {
//    using ScalarType = typename std::remove_cvref_t<decltype(uFunction)>::ctype;

//    using namespace Dune::DerivativeDirections;
//    using namespace Dune;
//
//    membraneStrainsAtVertices = std::vector<Eigen::Vector<ScalarType, 3>>();
//    for (int i = 0; auto& lP : lagrangePoints) {
//      const auto J                                = toEigen(geo.jacobianTransposed(lP));
//      const Eigen::Matrix<double, 2, 2> A         = J * J.transpose();
//      const Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(
//        uFunction.evaluateDerivative(lP, wrt(spatialAll, Dune::on(DerivativeDirections::referenceElement))));
//
//      const auto  epsV = defaultMembraneStrain.value(lP,geo,uFunction);
//
//      std::visit([&](auto& vec){
//        if constexpr (std::is_same_v<typename std::remove_cvref_t<decltype(vec)>::value_type,ScalarType>) {
//          vec.push_back(epsV);}
//          },membraneStrainsAtVertices);
//    }

  }
  template< typename Geometry>
  auto value(const Dune::FieldVector<double, 2> &gpPos,
             const Geometry &geo,
             const auto &uFunction) const -> Eigen::Vector3<typename std::remove_cvref_t<decltype(uFunction)>::ctype> {
    using ScalarType = typename std::remove_cvref_t<decltype(uFunction)>::ctype;

    cASAnsatzFunction.evaluateFunction(gpPos, NANS);

    Eigen::Vector<ScalarType, 3> res;
    res.setZero();
    for (int i = 0; auto& lP : lagrangePoints) {
      const auto  epsV = defaultMembraneStrain.value(lP,geo,uFunction);
      res += epsV * NANS[i++][0];
    }
//    std::visit([&](auto& vec){
//      if constexpr (std::is_same_v<typename std::remove_cvref_t<decltype(vec)>::value_type,ScalarType>) {
//        for (int i = 0; i < NANS.size(); ++i) {
//          res += vec[i] * NANS[i][0];
//        }
//      }
//      ;},membraneStrainsAtVertices);

    return res;
  }
  mutable Eigen::Matrix<double,Eigen::Dynamic,2> dN;
  template<typename Geometry,typename ScalarType>
  auto derivative(const Dune::FieldVector<double, 2> &gpPos,const Eigen::Matrix<ScalarType, 2, 3> &, const auto &, const Geometry& geo,const auto& uFunction, const auto& localBasis,
                  const int node) const{
    cASAnsatzFunction.evaluateFunction(gpPos, NANS);
    Eigen::Matrix<ScalarType, 3, 3> bop;
    bop.setZero();
    using namespace Dune::DerivativeDirections;
    using namespace Dune;

    for (int i = 0; auto& lP : lagrangePoints) {
       localBasis.evaluateJacobian(lP,dN);
//       std::cout<<"lP: "<<lP<<std::endl;
//       std::cout<<dN<<std::endl;


      const auto J = toEigen(geo.jacobianTransposed(lP));
      const Eigen::Matrix<ScalarType, 3, 2> gradu = toEigen(
        uFunction.evaluateDerivative(lP, wrt(spatialAll), Dune::on(DerivativeDirections::referenceElement)));
      const Eigen::Matrix<ScalarType, 2, 3> jE = J + gradu.transpose();
      const auto  bopI = defaultMembraneStrain.derivative<Geometry,ScalarType>(lP,jE,dN,geo,uFunction,localBasis,node);
      bop+=bopI*NANS[i++][0];
    }
    return bop;
  }

  template<typename Geometry,typename ScalarType>
  auto secondDerivative(const Dune::FieldVector<double, 2> &gpPos,const auto &, const Geometry& geo,const auto& uFunction, const auto& localBasis,
                        const Eigen::Vector3<ScalarType> &S, int I, int J)const {
    cASAnsatzFunction.evaluateFunction(gpPos, NANS);
    Eigen::Matrix<ScalarType, 3, 3> kg;
    kg.setZero();
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    for (int i = 0; auto& lP : lagrangePoints) {
     localBasis.evaluateJacobian(lP,dN);

      const auto  kgIJ = defaultMembraneStrain.secondDerivative(lP,dN,geo,uFunction,localBasis,S,I,J);
      kg+=kgIJ*NANS[i++][0];
    }
    return kg;
  }

};


template<class... Implementations>
class MembraneStrainVariant
{

 public:

  template<class Implementation>
  explicit MembraneStrainVariant(const Implementation& impl) :
      impl_(impl)
  {}

  MembraneStrainVariant() = default;
  MembraneStrainVariant(const MembraneStrainVariant& other) = default;
  template<class Implementation> requires (!std::is_same_v<Implementation,MembraneStrainVariant>)
  MembraneStrainVariant& operator=(const Implementation& impl)
      {
        impl_=impl;
        return *this;
      };
  MembraneStrainVariant(MembraneStrainVariant&& other)  noexcept = default;
  MembraneStrainVariant& operator=(const MembraneStrainVariant& other) = default;
  MembraneStrainVariant& operator=(MembraneStrainVariant&& other)  noexcept = default;

  template<typename Geometry>
  void pre( const Geometry &geo,
            const auto &uFunction){
  std::visit([&]( auto& impl) { impl.pre(geo,uFunction); }, impl_);
  }
  template< typename Geometry>
  auto value(const Dune::FieldVector<double, 2> &gpPos,
             const Geometry &geo,
             const auto &uFunction) {

    return std::visit([&](const auto& impl) { return impl.value(gpPos,geo,uFunction); }, impl_);
  }
  template<typename Geometry,typename ScalarType>
  auto derivative(const Dune::FieldVector<double, 2> &gpPos,const Eigen::Matrix<ScalarType, 2, 3> &jcur, const auto &dNAtGp, const Geometry& geo,const auto& uFunction, const auto& localBasis,
                  const int node) {

    return std::visit([&](const auto& impl) { return impl.derivative(gpPos,jcur,dNAtGp,geo,uFunction,localBasis,node); }, impl_);
  }
  template<typename Geometry,typename ScalarType>
  auto secondDerivative(const Dune::FieldVector<double, 2> &gpPos,const auto &dNAtGp, const Geometry& geo,const auto& uFunction, const auto& localBasis,
                        const Eigen::Vector3<ScalarType> &S, int I, int J) {

    return std::visit([&](const auto& impl) { return impl.secondDerivative(gpPos,dNAtGp,geo,uFunction,localBasis,S,I,J); }, impl_);
  }

 private:
  std::variant<Implementations...> impl_;
};

using MembraneStrain = MembraneStrainVariant<DefaultMembraneStrain,CASMembraneStrain<CASAnsatzFunction>,CASMembraneStrain<CASAnsatzFunctionANS>>;

}