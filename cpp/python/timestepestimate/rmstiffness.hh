// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "externalload.hh"

#include <utility>

//#include <dune/fufem/boundarypatch.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>
#include <dune/localfefunctions/cachedlocalBasis/cachedlocalBasis.hh>
#include <dune/localfefunctions/impl/standardLocalFunction.hh>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/fetraits.hh>
#include <ikarus/utils/eigendunetransformations.hh>
#include <ikarus/utils/linearalgebrahelper.hh>
#include <dune/geometry/quadraturerules.hh>
namespace Ikarus
{
  template <class BaseQuadrature, class Quadrature>
  auto tensorProductQuadrature(const BaseQuadrature & baseQuad, const Quadrature & onedQuad)
  {
    constexpr int baseQuadDim= BaseQuadrature::d;
    auto rule = Dune::QuadratureRule<double,baseQuadDim+1>() ;
    const unsigned int baseQuadSize = baseQuad.size();
    for( unsigned int bqi = 0; bqi < baseQuadSize; ++bqi )
    {
      const typename Dune::QuadraturePoint<double, baseQuadDim>::Vector &   basePoint = baseQuad[bqi].position( );
      const double &baseWeight = baseQuad[bqi].weight( );

      typename Dune::QuadraturePoint<double, baseQuadDim+1>::Vector point;
      for( unsigned int i = 0; i < baseQuadDim; ++i )
        point[ i ] = basePoint[ i ];

      const unsigned int onedQuadSize = onedQuad.size();
      for( unsigned int oqi = 0; oqi < onedQuadSize; ++oqi )
      {
        point[ baseQuadDim ] = onedQuad[oqi].position()[ 0 ];
        rule.emplace_back( Dune::QuadraturePoint(point, baseWeight * onedQuad[oqi].weight()) );
      }
    }
    return rule;
  }


  auto transformFromZeroOneToMinusOneToOne(const Dune::QuadraturePoint<double,1>& gp)
  {
    typename Dune::QuadraturePoint<double,1>::Vector mappedCoordinate = gp.position();
    mappedCoordinate[0]= 2*mappedCoordinate[0]-1.0;
    double mappedWeight = gp.weight() * 2.0;
    return Dune::QuadraturePoint<double,1>(mappedCoordinate, mappedWeight);
  }
}


namespace Ikarus {

  constexpr bool secondOrderBending = true;
  
  template<typename ResultantBasedShell>
  struct StressBasedShellRM : private ResultantBasedShell
  {
    template <typename ScalarType>
    using KinematicVariables =typename ResultantBasedShell::template KinematicVariables<ScalarType>;
    Dune::QuadratureRule<double,3> rule;
    static constexpr int midSurfaceDim         = 3;
    static constexpr int directorDim           = 3;
    static constexpr int directorCorrectionDim = directorDim - 1;
    template <typename... Args>
    explicit StressBasedShellRM(Args&&... args) : ResultantBasedShell{std::forward<Args>(args)...},geo{*this->geo_} {

      const auto& twoDRule = Dune::QuadratureRules<double,2>::rule(Dune::GeometryTypes::quadrilateral,this->order);
      const auto& oneDRule = Dune::QuadratureRules<double,1>::rule(Dune::GeometryTypes::line,4);
//      std::cout<<"oneDRule size "<<oneDRule.size()<<std::endl;
//      std::cout<<"oneDRule size "<<Dune::QuadratureRules<double,1>::rule(Dune::GeometryTypes::line,3).size()<<std::endl;
      numberOfThicknessIntegrationPoints = oneDRule.size();
      rule = Ikarus::tensorProductQuadrature(twoDRule,oneDRule);
//      fESettings = ResultantBasedShell::getfFESettings();
//      order = ResultantBasedShell::getOrder();
//      numberOfNodes = ResultantBasedShell::getNumberOfNodes();
//      membraneStrain = ResultantBasedShell::getMembraneStrain();
//      std::cout<<"Constructor: "<<rule.size()<<std::endl;
      secondOrderBending = true;

    }

    struct FirstOrderVariationVariables
    {

    };


    int numberOfThicknessIntegrationPoints;
    bool secondOrderBending=true;
//    FESettings fESettings;
//    size_t numberOfNodes{0};
//    int order{};
//    mutable MembraneStrain membraneStrain;

    using GlobalIndex = typename ResultantBasedShell::GlobalIndex;
    using ResultantBasedShell::size ;
    using ResultantBasedShell::globalFlatIndices ;
    using ResultantBasedShell::localView ;
    using ResultantBasedShell::ndofEmbedded ;
    using Basis = typename ResultantBasedShell::Basis;
    using FlatBasis = typename Basis::FlatBasis;
    using FERequirementType = typename ResultantBasedShell::FERequirementType;
    using ResultRequirementsType = ResultRequirements<FERequirementType>;
    using LocalView = typename FlatBasis::LocalView;
    using Element = typename LocalView::Element;
    using Geometry = typename Element::Geometry;
    const Geometry& geo;
    using GridView = typename FlatBasis::GridView;
    static constexpr int useEigenRef = ResultantBasedShell::useEigenRef;
    using Traits = FETraits<Element,useEigenRef>;
   protected:
    static constexpr int myDim = Traits::mydim;
    static constexpr int worlddim = Traits::worlddim;
    static constexpr bool isEuclidean = ResultantBasedShell::isEuclidean;
    using ResultantBasedShell::localViewBlocked ;
   public:
    inline void calculateMatrix(const FERequirementType &par, typename Traits::MatrixType K) const {
      calculateMatrixImpl<double>(par, K);
    }

    inline void calculateVector(const FERequirementType &par, typename Traits::VectorType force) const {
      calculateVectorImpl<double>(par, force);
    }
    inline double calculateScalar(const FERequirementType &par) const { return calculateScalarImpl<double>(par); }




    template<typename ScalarType>
    auto calculateScalarImpl(const FERequirementType &par, const std::optional<const Eigen::VectorX<ScalarType>> &dx
                                                           = std::nullopt) const -> ScalarType {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto [ displacementFunction, directorFunction,
          directorReferenceFunction]
          = this->createFunctions(par,dx);
      const auto &lambda = par.getParameter(FEParameter::loadfactor);
      ScalarType energy = 0.0;
      this->membraneStrain.pre(geo,displacementFunction);

      const auto &thickness_ = this->fESettings.template request<double>("thickness");
      KinematicVariables<ScalarType> kin{};
      for (int gpIndex=0; const auto & gp: rule) {
        const double zeta  = (2*gp.position()[2]-1)*thickness_/2.0;
        const int gpIndex2D= gpIndex/numberOfThicknessIntegrationPoints;
        const Dune::FieldVector<double,2> gp2DPos= {gp.position()[0],gp.position()[1]};

        kin.t           = directorFunction.evaluate(gpIndex2D,Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.t0          = directorReferenceFunction.evaluate(gpIndex2D,Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.ud1andud2   = displacementFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.A1andA2     = Dune::toEigen(geo.jacobianTransposed(gp2DPos)).transpose();
        kin.a1anda2     = kin.A1andA2+ kin.ud1andud2;
        kin.t0d1Andt0d2 = directorReferenceFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.td1Andtd2   = directorFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));

        auto [_,epsV,kappaV,gammaV]                = this->computeMaterialAndStrains(gp2DPos,gpIndex2D,geo,displacementFunction,directorFunction,directorReferenceFunction,kin);

//        const auto
//            [C, epsV, kappaV, j, J, h,H, a3N, a3] = this->computeMaterialAndStrains(gp2DPos, geo, uFunction);

        const auto G = this->calc3DMetric(kin,zeta);

        const auto intElement =sqrt(G.determinant()) * gp.weight() * 2 * thickness_ / 2.0;

//        std::cout<<"G: "<<G<<std::endl;
        const auto Ginv = G.inverse().eval();

        const auto C3D = this->materialTangent(Ginv);
        Eigen::Vector3<ScalarType> rhoV;

        rhoV<< 0.5*(kin.td1().squaredNorm()-kin.t0d1().squaredNorm()),0.5*(kin.td2().squaredNorm()-kin.t0d2().squaredNorm()),kin.td1().dot(kin.td2())-kin.t0d1().dot(kin.t0d2());
        const auto strainsV= (epsV+ zeta*kappaV+secondOrderBending*zeta*zeta*rhoV).eval();
        Eigen::Vector<ScalarType,5> strains;
        strains<< strainsV,gammaV;
        const ScalarType energyVal = 0.5*strains.dot(C3D*strains);
        energy += (energyVal)*intElement;
        ++gpIndex;
      }

      if (this->volumeLoad) {
        for (const auto& [gpIndex, gp] : displacementFunction.viewOverIntegrationPoints()) {
          const auto u                                       = displacementFunction.evaluate(gpIndex);
          const Eigen::Vector<double, Traits::worlddim> fExt = this->volumeLoad(geo.global(gp.position()), lambda)[0];
          energy -= u.dot(fExt) * geo.integrationElement(gp.position()) * gp.weight();
        }
      }

      // line or surface loads, i.e., neumann boundary
      if (not this->neumannBoundary and not this->neumannBoundaryLoad) return energy;
      forEachInterSectionIntegrationPoint(this->localView().element(),this->neumannBoundary,this->order,
                                          [&](auto& quadPos,auto&& globalPos,auto& intElement){
                                            const auto [fext,mext]
                                                = this->neumannBoundaryLoad(globalPos, lambda);
                                            const auto u = displacementFunction.evaluate(quadPos);
//                                            const auto t = directorFunction.evaluate(quadPos);
                                            energy -= fext.dot(u) * intElement;
//                                            energy -= mext.dot(t) * intElement;
                                          });
      return energy;
    }

    template<typename ScalarType>
    void calculateVectorImpl(const FERequirementType &par, typename Traits::VectorType force,
                             const std::optional<const Eigen::VectorX<ScalarType>> &dx = std::nullopt) const {
      force.setZero();
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto [ displacementFunction, directorFunction,
          directorReferenceFunction]
          = this->createFunctions(par,dx);
      const auto &lambda = par.getParameter(FEParameter::loadfactor);
      const int midSurfaceDofs = this->numNodeMidSurface * midSurfaceDim;
      const auto &thickness_ = this->fESettings.template request<double>("thickness");
      this->membraneStrain.pre(geo,displacementFunction);
      KinematicVariables<ScalarType> kin{};
      // Internal forces
      Eigen::Matrix<ScalarType, 5, 3> bopMidSurfaceI;
      Eigen::Matrix<ScalarType, 5, 2> bopDirectorI;

      for (int gpIndex=0; const auto & gp: rule) {
        const double zeta  = (2*gp.position()[2]-1)*thickness_/2.0;
        const int gpIndex2D= gpIndex/numberOfThicknessIntegrationPoints;

        const Dune::FieldVector<double,2> gp2DPos= {gp.position()[0],gp.position()[1]};
        kin.t           = directorFunction.evaluate(gpIndex2D,Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.t0          = directorReferenceFunction.evaluate(gpIndex2D,Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.ud1andud2   = displacementFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.A1andA2     = Dune::toEigen(geo.jacobianTransposed(gp2DPos)).transpose();
        kin.a1anda2     = kin.A1andA2+ kin.ud1andud2;
        kin.t0d1Andt0d2 = directorReferenceFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.td1Andtd2   = directorFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        const Eigen::Matrix<ScalarType,2,3> jE = kin.a1anda2.transpose();
        auto [_,epsV,kappaV,gammaV]                = this->computeMaterialAndStrains(gp2DPos,gpIndex2D,geo,displacementFunction,directorFunction,directorReferenceFunction,kin);

        const auto G = this->calc3DMetric(kin,zeta);
        //        std::cout<<"G: "<<G<<std::endl;
        const auto Ginv = G.inverse().eval();

        const auto g1Andg2 = (kin.a1anda2+ zeta* kin.td1Andtd2).eval();

        const auto intElement =sqrt(G.determinant()) * gp.weight() * 2 * thickness_ / 2.0;


        const auto C3D = this->materialTangent(Ginv);
        Eigen::Vector3<ScalarType> rhoV;
        rhoV<< 0.5*(kin.td1().squaredNorm()-kin.t0d1().squaredNorm()),0.5*(kin.td2().squaredNorm()-kin.t0d2().squaredNorm()),kin.td1().dot(kin.td2())-kin.t0d1().dot(kin.t0d2());
        const auto strainsV= (epsV+ zeta*kappaV+secondOrderBending*zeta*zeta*rhoV).eval();
        Eigen::Vector<ScalarType,5> strains;
        strains<< strainsV,gammaV;
        const auto S = (C3D*strains).eval();

        const auto &Nd = this->localBasisMidSurface.evaluateJacobian(gpIndex2D);
        for (size_t i = 0; i < this->numNodeMidSurface; ++i) {
          bopMidSurfaceI= bopMembrane(kin, gpIndex2D, gp2DPos,jE, Nd, displacementFunction,directorFunction, i, zeta);

          force.template segment<3>(3*i) +=
              bopMidSurfaceI.transpose()*S*intElement;
          // the first two fixes the change of the integration mapping from 0..1 to -1..1,
          // and the h/2 factor is the factor for the correct thickness
        }
        for (int i = 0; i < this->numNodeDirector; ++i) {
          const auto indexI = midSurfaceDofs + directorCorrectionDim * i;

          bopDirectorI =secondOrderBending? bopDirector(kin,Nd, g1Andg2, gpIndex2D, gp2DPos, displacementFunction, directorFunction, i, zeta)
                                            :bopDirector(kin,Nd, kin.a1anda2, gpIndex2D, gp2DPos,displacementFunction, directorFunction, i, zeta);
          force.template segment<directorCorrectionDim>(indexI) += bopDirectorI.transpose() * S*intElement;
        }

        ++gpIndex;
      }
      using namespace Dune::Indices;
      //External forces volume forces over the domain
      if (this->volumeLoad) {
        for (const auto& [gpIndex, gp] : displacementFunction.viewOverIntegrationPoints()) {
          const auto [fext,mext] = this->volumeLoad(geo.global(gp.position()), lambda);
          for (size_t i = 0; i < this->numNodeMidSurface; ++i) {
            const auto udCi = displacementFunction.evaluateDerivative(gpIndex, wrt(coeff(i)));
            force.template segment<worlddim>(worlddim * i)
                -= udCi * fext * geo.integrationElement(gp.position()) * gp.weight();
          }
//          for (size_t i = 0; i < this->numNodeDirector; ++i) {
//            const auto indexI = midSurfaceDofs + directorCorrectionDim * i;
//            const auto tdCi = directorFunction.evaluateDerivative(gpIndex, Dune::wrt(coeff(_1,i)));
//            force.template segment<directorCorrectionDim>(indexI)
//                -= (tdCi.transpose() * mext) * geo.integrationElement(gp.position()) * gp.weight();
//          }
        }
      }

      // External forces, boundary forces, i.e., at the Neumann boundary
      if (not this->neumannBoundary) return;

      forEachInterSectionIntegrationPoint(this->localView().element(),this->neumannBoundary,this->order,
                                          [&](auto& quadPos,auto&& globalPos,auto& intElement){
                                            const auto [fext,mext] = this->neumannBoundaryLoad(globalPos, lambda);

                                            for (size_t i = 0; i < this->numNodeMidSurface; ++i) {
                                              const auto udCi = displacementFunction.evaluateDerivative(quadPos, wrt(coeff(i)));

                                              // Value of the Neumann data at the current position
                                              force.template segment<worlddim>(worlddim * i) -= udCi * fext * intElement;
                                            }
//                                            const auto &Ndirector = this->localBasisDirector.evaluateFunction(quadPos);
                                            for (size_t i = 0; i < this->numNodeDirector; ++i) {
                                              const auto indexI = midSurfaceDofs + directorCorrectionDim * i;
                                              const auto tdCi = directorFunction.evaluateDerivative(quadPos, Dune::wrt(coeff(_1,i)));
                                              const auto t=  directorFunction.evaluate(quadPos,Dune::on(Dune::DerivativeDirections::referenceElement));
                                              force.template segment<directorCorrectionDim>(indexI)
                                                  -= (tdCi.transpose() * mext.cross(t)) * intElement;
                                            }
                                          });
    }
    template<typename ScalarType>
    auto bopMembrane(const KinematicVariables<ScalarType>& kin, const auto& gpIndex2D,const auto& gp2DPos,const auto& jE, const auto& Nd, const auto& displacementFunction, const auto& directorFunction, int i,double zeta) const
    {
      Eigen::Matrix<ScalarType, 5, 3> bopMidSurfaceI;
      const auto bopIMembrane   = this->membraneStrain.derivative(gp2DPos,jE, Nd, geo,displacementFunction,this->localBasisMidSurface, i);
      const auto bopIBending   = this->boperatorMidSurfaceBending(kin,gpIndex2D, i,displacementFunction);
//      const auto bopIShear   = this->boperatorMidSurfaceShear(kin,gpIndex2D, i,displacementFunction);
      const auto bopIShear   = this-> transverseShearStrain.derivativeWRTMidSurface(kin,gp2DPos,gpIndex2D, Nd,geo,displacementFunction,directorFunction,this->localBasisMidSurface,i);
      bopMidSurfaceI<< bopIMembrane+zeta*bopIBending, bopIShear;
      return bopMidSurfaceI;
    }
    template<typename ScalarType>
    auto bopDirector(const KinematicVariables<ScalarType>& kin,const auto& Nd, const auto& g1Andg2,
                     const auto& gpIndex2D,const auto& gp2DPos,const auto& displacementFunction,
                     const auto& directorFunction,int i,double zeta)const
    {
      Eigen::Matrix<ScalarType, 5, 2> bopDirectorI;

      const auto bopBendingI   = this->boperatorDirectorBending(g1Andg2, gpIndex2D, i, directorFunction);
      const auto bopShearI   = this->transverseShearStrain.derivativeWRTDirector(kin,gp2DPos,gpIndex2D, Nd,geo,displacementFunction,directorFunction,
                                                                                   this->localBasisMidSurface,i);
      bopDirectorI<<zeta*bopBendingI,bopShearI;
      return bopDirectorI;

    }
    template<typename ScalarType>
    void calculateMatrixImpl(const FERequirementType &par, typename Traits::MatrixType K,
                             const std::optional<const Eigen::VectorX<ScalarType>> &dx = std::nullopt) const {
    
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto [ displacementFunction, directorFunction,
          directorReferenceFunction]
          = this->createFunctions(par,dx);
      const auto &lambda = par.getParameter(FEParameter::loadfactor);
      this->membraneStrain.pre(geo,displacementFunction);
      const int midSurfaceDofs = this->numNodeMidSurface * midSurfaceDim;
      KinematicVariables<ScalarType> kin{};

std::cout<<"Step1"<<std::endl;
      const auto &thickness_ = this->settings_.thickness;
      Eigen::Matrix<ScalarType, 5, 3> bopMidSurfaceI;
      Eigen::Matrix<ScalarType, 5, 2> bopDirectorI;
      Eigen::Matrix<ScalarType, 5, 3> bopMidSurfaceJ;
      Eigen::Matrix<ScalarType, 5, 2> bopDirectorJ;
      // Internal forces
      for (int gpIndex=0; const auto & gp: rule) {
        const double zeta  = (2*gp.position()[2]-1)*thickness_/2.0;
        const int gpIndex2D= gpIndex/numberOfThicknessIntegrationPoints;
        const Dune::FieldVector<double,2> gp2DPos= {gp.position()[0],gp.position()[1]};
        kin.t           = directorFunction.evaluate(gpIndex2D,Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.t0          = directorReferenceFunction.evaluate(gpIndex2D,Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.ud1andud2   = displacementFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.A1andA2     = Dune::toEigen(geo.jacobianTransposed(gp2DPos)).transpose();
        kin.a1anda2     = kin.A1andA2+ kin.ud1andud2;
        kin.t0d1Andt0d2 = directorReferenceFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.td1Andtd2   = directorFunction.evaluateDerivative(gpIndex2D, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        const Eigen::Matrix<ScalarType,2,3> jE = kin.a1anda2.transpose();
        auto [_,epsV,kappaV,gammaV]                = this->computeMaterialAndStrains(gp2DPos,gpIndex2D,geo,displacementFunction,directorFunction,directorReferenceFunction,kin);

        const auto G = this->calc3DMetric(kin,zeta);
        //        std::cout<<"G: "<<G<<std::endl;
        const auto intElement =sqrt(G.determinant()) * gp.weight() * 2 * thickness_ / 2.0;

        const auto Ginv = G.inverse().eval();

        const auto g1Andg2 = (kin.a1anda2+ zeta* kin.td1Andtd2).eval();

        const auto C3D = this->materialTangent(Ginv);
        Eigen::Vector3<ScalarType> rhoV;
        rhoV<< 0.5*(kin.td1().squaredNorm()-kin.t0d1().squaredNorm()),0.5*(kin.td2().squaredNorm()-kin.t0d2().squaredNorm()),kin.td1().dot(kin.td2())-kin.t0d1().dot(kin.t0d2());
        const auto strainsV= (epsV+ zeta*kappaV+secondOrderBending*zeta*zeta*rhoV).eval();
        Eigen::Vector<ScalarType,5> strains;
        strains<< strainsV,gammaV;
        const auto S = (C3D*strains).eval();
        Eigen::Vector<ScalarType,8> S8;
        Eigen::Vector<ScalarType,3> SSec;
        SSec= zeta*zeta*S.template segment<3>(0);
        S8.template segment<3>(0)<< S.template segment<3>(0);
        S8.template segment<3>(3)<< zeta*S.template segment<3>(0);
        S8.template segment<2>(6)<< S.template segment<2>(3);
        const auto &Nd = this->localBasisMidSurface.evaluateJacobian(gpIndex2D);
        const auto &N = this->localBasisMidSurface.evaluateFunction(gpIndex2D);
        const auto &dNdirector = this->localBasisDirector.evaluateJacobian(gpIndex2D);
        const auto &Ndirector = this->localBasisDirector.evaluateFunction(gpIndex2D);
        for (size_t i = 0; i < this->numNodeMidSurface; ++i) {
          const auto indexI = midSurfaceDim * i;
          bopMidSurfaceI= bopMembrane(kin, gpIndex2D, gp2DPos,jE, Nd, displacementFunction,directorFunction, i, zeta);

          for (size_t j = i; j < this->numNodeMidSurface; ++j) {
            const auto indexJ = midSurfaceDim * j;

            bopMidSurfaceJ= bopMembrane(kin, gpIndex2D, gp2DPos,jE, Nd, displacementFunction,directorFunction, j, zeta);

            K.template block<3, 3>(indexI, indexJ) += (bopMidSurfaceI.transpose()*C3D*bopMidSurfaceJ)*intElement;

            Eigen::Matrix<ScalarType, 3, 3> kgMembraneIJ = this->membraneStrain.secondDerivative(gp2DPos,Nd,geo,displacementFunction,this->localBasisMidSurface, S.template segment<3>(0).eval(), i, j);
            K.template block<3, 3>(indexI, indexJ) += kgMembraneIJ*intElement;

          }
          for (int j = 0; j < this->numNodeDirector; ++j) {
            const auto indexJ = midSurfaceDofs + directorCorrectionDim * j;
//            const auto bopBendingJ   = this->boperatorDirectorBending(g1Andg2, gpIndex2D, j, directorFunction);
//            const auto bopShearJ   = this->boperatorDirectorShear(kin, gpIndex2D, j, directorFunction);
            bopDirectorJ =secondOrderBending? bopDirector(kin,Nd, g1Andg2, gpIndex2D,gp2DPos, displacementFunction, directorFunction, j, zeta):bopDirector(kin,Nd, kin.a1anda2, gpIndex2D,gp2DPos, displacementFunction,directorFunction, j, zeta);

            Eigen::Matrix<ScalarType, 3, 2> kg= this->kgMidSurfaceDirectorBending(kin,N,Nd,gpIndex2D,i,j,displacementFunction,directorFunction,S8);
//            Eigen::Matrix<ScalarType, 3, 2> kg2= this->kgMidSurfaceDirectorShear(kin,N,Nd,gpIndex2D,i,j,displacementFunction,directorFunction,S8);
            Eigen::Matrix<ScalarType, 3, 2> kg2= this->transverseShearStrain.secondDerivativeWRTSurfaceDirector(gp2DPos,gpIndex, Nd,geo,displacementFunction,directorFunction,
                                                                                                                 this->localBasisMidSurface,
                                                                                                                 S8,kin,i,j);

            K.template block<midSurfaceDim, directorCorrectionDim>(indexI, indexJ)
                += bopMidSurfaceI.transpose() * C3D * bopDirectorJ * intElement;

            K.template block<midSurfaceDim, directorCorrectionDim>(indexI, indexJ)
                += (kg+kg2) * intElement;
          }
        }
        for (int i = 0; i < this->numNodeDirector; ++i) {
          const auto indexI = midSurfaceDofs + directorCorrectionDim * i;
          bopDirectorI =secondOrderBending? bopDirector(kin,Nd, g1Andg2, gpIndex2D,gp2DPos,displacementFunction, directorFunction, i, zeta) : bopDirector(kin,Nd, kin.a1anda2, gpIndex2D,gp2DPos,displacementFunction, directorFunction, i, zeta);

          for (int j = i; j < this->numNodeDirector; ++j) {
            const auto indexJ = midSurfaceDofs + directorCorrectionDim * j;

            bopDirectorJ =secondOrderBending? bopDirector(kin,Nd, g1Andg2, gpIndex2D,gp2DPos,displacementFunction, directorFunction, j, zeta): bopDirector(kin,Nd, kin.a1anda2, gpIndex2D,gp2DPos,displacementFunction, directorFunction, j, zeta);

            Eigen::Matrix<ScalarType, 2, 2> kgBending= this->kgDirectorDirectorBending(kin,Ndirector,dNdirector,gpIndex2D,i,j,displacementFunction,directorFunction,S8);
            Eigen::Matrix<ScalarType, 2, 2> kgBending2= this->kgSecondDirectorDirectorBending(kin,Ndirector,dNdirector,gpIndex2D,i,j,displacementFunction,directorFunction,SSec);
//            Eigen::Matrix<ScalarType, 2, 2> kgShear= this->kgDirectorDirectorShear(kin,Ndirector,dNdirector,gpIndex2D,i,j,displacementFunction,directorFunction,S8);
            Eigen::Matrix<ScalarType, 2, 2> kgShear= this->transverseShearStrain.secondDerivativeWRTDirectorDirector(gp2DPos,gpIndex, Nd,geo,displacementFunction,directorFunction,
                                                                                                                        this->localBasisMidSurface,S8,kin,i,j);

            K.template block<directorCorrectionDim, directorCorrectionDim>(indexI, indexJ)
                += bopDirectorI.transpose() * C3D * bopDirectorJ * intElement;

            K.template block<directorCorrectionDim, directorCorrectionDim>(indexI, indexJ)
                += (kgBending+secondOrderBending*kgBending2+kgShear)* intElement;
          }
        }
        ++gpIndex;
      }

      K.template triangularView<Eigen::StrictlyLower>() = K.transpose();
          }

    auto calculateStresses(const auto& par,
        const Dune::FieldVector<double,3>& gp, bool toIntegrate
    ) const
    {
//       std::array<Eigen::Matrix<double,3,3>,3> stresses; //Cauchy,PK2,PK1

//        KinematicVariables<double> kin;
//       const auto &thickness_ = this->fESettings.template request<double>("thickness");
//       using namespace Dune::DerivativeDirections;
//       using namespace Dune;
//       const auto [ displacementFunction, directorFunction,
//           directorReferenceFunction]
//           = this->template createFunctions<double>(par);
//       const auto &lambda = par.getParameter(FEParameter::loadfactor);
//       const double zeta  = (2*gp[2]-1)*thickness_/2.0;
//       const Dune::FieldVector<double,2> gp2DPos= {gp[0],gp[1]};
//       kin.t           = directorFunction.evaluate(gp2DPos,Dune::on(Dune::DerivativeDirections::referenceElement));
//       kin.t0          = directorReferenceFunction.evaluate(gp2DPos,Dune::on(Dune::DerivativeDirections::referenceElement));
//       kin.ud1andud2   = displacementFunction.evaluateDerivative(gp2DPos, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
//       kin.A1andA2     = Dune::toEigen(geo.jacobianTransposed(gp2DPos)).transpose();
// //      kin.A1andA2     = Eigen::Matrix<double,3,2>::Identity();
// //      kin.a1anda2     = kin.A1andA2+ kin.ud1andud2;
//       kin.t0d1Andt0d2 = directorReferenceFunction.evaluateDerivative(gp2DPos, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
//       kin.td1Andtd2   = directorFunction.evaluateDerivative(gp2DPos, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
//       using namespace Dune::Indices;
//       const auto &mNodal = par.getGlobalSolution(Ikarus::FESolutions::noSolution)[_0];
//     //      const auto &lambda = par.getParameter(Ikarus::FEParameter::loadfactor);

//       const auto &child0 = this->localViewBlocked_.tree().child(_0, 0);
//       const auto &fe0    = child0.finiteElement();

//       Dune::BlockVector<typename std::remove_cvref_t<decltype(mNodal[0])>>  displacements(fe0.size());
//       const auto& cps= geo.impl().controlPoints().directGetAll();
//       for (auto i = 0U; i < fe0.size(); ++i) {
//         const auto globalIndex = this->localViewBlocked_.index(this->localViewBlocked_.tree().child(_0).localIndex(i));
//         displacements[i].setValue(mNodal[globalIndex[1]].getValue() + toEigen(cps[i]));
//       }
// //      std::cout<<displacements<<std::endl;
//       Dune::StandardLocalFunction midSurface(this->localBasisMidSurface, displacements, this->geo_, _0);
//       Eigen::MatrixX2d Nd;
//       this->localBasisMidSurface.evaluateJacobian(gp2DPos,Nd);

//       const auto a1anda2T= (Dune::viewAsEigenMatrixAsDynFixed(displacements).transpose()*Nd).eval();

//       kin.a1anda2   = midSurface.evaluateDerivative(gp2DPos, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));

//       const Eigen::Matrix<double,2,3> jE = kin.a1anda2.transpose();
//       const int dummyIndex=-100;
//       auto [_,epsV,kappaV,gammaV]                = this->computeMaterialAndStrains(gp2DPos,dummyIndex,geo,displacementFunction,directorFunction,directorReferenceFunction,kin);
// //      std::cout<<"gammaVInside"<<std::endl;
// //      std::cout<<gammaV<<std::endl;
//       const auto G = this->calc3DMetric(kin,zeta);
//       const auto Ginv = G.inverse().eval();

//       const auto g1Andg2 = (kin.a1anda2+ zeta* kin.td1Andtd2).eval();
//       const auto G1AndG2 = (kin.A1andA2+ zeta* kin.t0d1Andt0d2).eval();

//       const auto C3D = this->materialTangent(Ginv);
//       Eigen::Vector3<double> rhoV;
//       rhoV<< 0.5*(kin.td1().squaredNorm()-kin.t0d1().squaredNorm()),0.5*(kin.td2().squaredNorm()-kin.t0d2().squaredNorm()),kin.td1().dot(kin.td2())-kin.t0d1().dot(kin.t0d2());
//       const auto strainsV= (epsV+ zeta*kappaV+secondOrderBending*zeta*zeta*rhoV).eval();
//       Eigen::Vector<double,5> strains;
//       Eigen::Vector<double,6> strains6;
//       strains<< strainsV,gammaV;
//       strains6<<strains,0.0;
// //      std::cout<<"strains6"<<std::endl;
// //      std::cout<<strains6<<std::endl;

//       Eigen::Matrix3d E;
//       E<< strains(0),strains(2)/2, strains(3)/2,
//           strains(2)/2, strains(1), strains(4)/2,
//           strains(3)/2,strains(4)/2,0;

//       Eigen::Matrix3d PK2 ;

//       Eigen::Matrix3d J;
//       Eigen::Matrix3d j;
//       J<<G1AndG2.col(0),G1AndG2.col(1),kin.t0;
//       j<<g1Andg2.col(0),g1Andg2.col(1),kin.t;

//       const  Eigen::Matrix3d Jloc= orthonormalizeMatrixColumns(J);
//       const  Eigen::Matrix3d jloc= orthonormalizeMatrixColumns(j);
//       const  Eigen::Matrix3d F     = j*(J.inverse());

//       const double detF = F.determinant();

//       const auto &emod_ = this->fESettings.template request<double>("youngs_modulus");
//       const auto &nu_ = this->fESettings.template request<double>("poissons_ratio");
//       const double lambdaf = emod_*nu_/((1.0 + nu_)*(1.0 - 2.0*nu_));
//       const double mu = emod_/(2.0*(1.0 + nu_));
//       const double lambdbar = 2.0*lambdaf*mu/(lambdaf + 2.0*mu);
// //       Eigen::Matrix3d E = 0.5*(F.transpose()*F-Eigen::Matrix3d::Identity());
// //       if(not secondOrderBending)
// //       E=Ikarus::fromVoigt(strains6,true);
//        const Eigen::Matrix3d T= Jloc*J.inverse();
// //       std::cout<<"EB"<<std::endl;
// //       std::cout<<E<<std::endl;
//        E= T.transpose()*E*T;
// //       std::cout<<"E"<<std::endl;
// //       std::cout<<E<<std::endl;
//       PK2= lambdbar* E.trace()*Eigen::Matrix3d::Identity()+2*mu*E;
// //      std::cout<<"PK2"<<std::endl;
// //      std::cout<<PK2<<std::endl;
//       Eigen::Matrix3d cauchy=1/detF*F*PK2*F.transpose();
// //      std::cout<<"cauchyB"<<std::endl;
// //      std::cout<<cauchy<<std::endl;
// //      if(toIntegrate)
// //      {
// //        Eigen::Matrix3d jC;
// //        jC<<      kin.a1anda2.col(0),kin.a1anda2.col(1),kin.t;
// ////        JC<< zujhjhujh    kin.A1andA2.col(0),kin.A1andA2.col(1),kin.t0;
// //
// //        const  Eigen::Matrix3d jlocC= orthonormalizeMatrixColumns(jC);
// //        cauchy=jlocC.transpose()*cauchy*jlocC;
// //
// //      }else {
//         cauchy = jloc.transpose() * cauchy * jloc;
// //      }
// //      std::cout<<"cauchyA"<<std::endl;
// //      std::cout<<cauchy<<std::endl;
//       stresses[0]= cauchy;
//       stresses[1]= PK2;
//       double E11 = E(0,0);
//       return std::make_tuple(stresses,detF,E11,zeta,j,J,gp2DPos);
    }

    void calculateAt(const ResultRequirementsType &req, const Dune::FieldVector<double, 3> &local,
                     ResultTypeMap<double> &result) const {
      //   typename ResultTypeMap<double>::ResultArray resultVector;

      //   resultVector.resize(3, 3);
      //   auto [res,detF,E11,zeta,j,J,gp2DPos]= calculateStresses(req.getFERequirements(),local,false);
      //   if(req.isResultRequested(ResultType::cauchyStress)) {
      //     resultVector = res[0];
      //     result.insertOrAssignResult(ResultType::cauchyStress, resultVector);
      //   }else if(req.isResultRequested(ResultType::PK2Stress)) {
      //     resultVector = res[1];
      //     result.insertOrAssignResult(ResultType::PK2Stress, resultVector);
      //   }else if(req.isResultRequested(ResultType::detF)) {
      //     resultVector.resize(1,1);
      //     resultVector(0,0) = detF;
      //     result.insertOrAssignResult(ResultType::detF, resultVector);
      //   }else if(req.isResultRequested(ResultType::E11)) {
      //     resultVector.resize(1,1);
      //     resultVector(0,0) = E11;
      //     result.insertOrAssignResult(ResultType::E11, resultVector);
      //   }else if(req.isResultRequested(ResultType::zeta)) {
      //     resultVector.resize(1,1);
      //     resultVector(0,0) = zeta;
      //     result.insertOrAssignResult(ResultType::zeta, resultVector);
      //   }else if(req.isResultRequested(ResultType::refJacobian)) {
      //     resultVector.resize(3,3);
      //     resultVector = J;
      //     result.insertOrAssignResult(ResultType::refJacobian, resultVector);
      //   }else if(req.isResultRequested(ResultType::curJacobian)) {
      //     resultVector.resize(3,3);
      //     resultVector = j;
      //     result.insertOrAssignResult(ResultType::curJacobian, resultVector);
      //   }else if(req.isResultRequested(ResultType::gpPosition)) {
      //     resultVector.resize(2,1);
      //     resultVector(0,0) = gp2DPos[0];
      //     resultVector(1,0) = gp2DPos[1];
      //     result.insertOrAssignResult(ResultType::gpPosition, resultVector);
      //   }
      // else
      // DUNE_THROW(Dune::NotImplemented, "No results are implemented");
    }

    void calculateAt(const ResultRequirementsType &req, const Dune::FieldVector<double, 2> &local,
                     ResultTypeMap<double> &result) const {
      // if(req.isResultRequested(ResultType::energyArray)) {
      //   typename ResultTypeMap<double>::ResultArray resultVector;
      //   auto energyVec = calculateInternalEnergy(req.getFERequirements(),local);
      //   resultVector.resize(3, 1);
      //   resultVector=energyVec;
      //   result.insertOrAssignResult(ResultType::energyArray, resultVector);
      //   return;
      // }
      // if (req.isResultRequested(ResultType::cauchyStress))
      //   DUNE_THROW(Dune::NotImplemented, "No results are implemented");
      // else if (req.isResultRequested(ResultType::membraneForces)
      //          or req.isResultRequested(ResultType::bendingMoments)
      //          or req.isResultRequested(ResultType::shearForces)or req.isResultRequested(ResultType::membraneForcesPK2)
      //          or req.isResultRequested(ResultType::bendingMomentsPK2)
      //          or req.isResultRequested(ResultType::shearForcesPK2))
      // {
      //    stressResultants(req,local,result);
      // }else
      //   DUNE_THROW(Dune::NotImplemented, "No results are implemented");
    }

    auto calculateEnergyPart(const ResultRequirementsType &req, const Dune::FieldVector<double, 2> &local,
                             ResultTypeMap<double> &result)
    {

    }

    auto calcStressResultants(const FERequirementType &par, const Dune::FieldVector<double, 2> &local)const
    {

      // const auto& oneDRule = Dune::QuadratureRules<double,1>::rule(Dune::GeometryTypes::line,4);
      // const auto &thickness = this->fESettings.template request<double>("thickness");
      // using namespace Dune::DerivativeDirections;
      // using namespace Dune;
      // const auto [ displacementFunction, directorFunction,
      //     directorReferenceFunction]
      //     = this->template createFunctions<double>(par);

      // Eigen::Matrix2d NCauchy = Eigen::Matrix2d::Zero();
      // Eigen::Matrix2d MCauchy = Eigen::Matrix2d::Zero();
      // Eigen::Matrix2d MCauchy2 = Eigen::Matrix2d::Zero();
      // Eigen::Vector2d QCauchy = Eigen::Vector2d::Zero();

      // Eigen::Matrix2d NPK2 = Eigen::Matrix2d::Zero();
      // Eigen::Matrix2d MPK2 = Eigen::Matrix2d::Zero();
      // Eigen::Matrix2d MPK2_2 = Eigen::Matrix2d::Zero();
      // Eigen::Vector2d QPK2 = Eigen::Vector2d::Zero();

      // KinematicVariables<double> kin;
      // for (auto& gP : oneDRule)
      // {
      //   const auto& gppos =gP.position();
      //   const auto& gpweight =gP.weight();

      //   const double zeta  = (2*gppos[0]-1)*thickness/2.0;
      //   const Dune::FieldVector<double,2> gp2DPos= local;
      //   const Dune::FieldVector<double,3> gp3DPos= {local[0],local[1],gppos[0]};

      //   const double fac_gp = gpweight * thickness * 0.5*2;
      //   // the first two fixes the change of the integration mapping from 0..1 to -1..1,
      //   // and the h/2 factor is the factor for the correct thickness

      //   const auto [res,detF,E11,zetaS,jS,JS,gp2DPosS]   = calculateStresses(par,gp3DPos,true);

      //   const auto [ cauchy,PK2,PK1] = res;
      //   Eigen::Matrix2d PK2_al_be= PK2.template block<2,2>(0,0);
      //   Eigen::Matrix2d cauchy_al_be= cauchy.template block<2,2>(0,0);

      //   const Eigen::Vector2d cauchy_al_3              = cauchy.template block<2,1>(0,2);
      //   const Eigen::Vector2d PK2_al_3              = PK2.template block<2,1>(0,2);

      //   NCauchy +=  fac_gp *  cauchy_al_be ;
      //   MCauchy +=  fac_gp *  zeta *cauchy_al_be ;
      //   MCauchy2 +=  fac_gp *  zeta * zeta *cauchy_al_be ;
      //   QCauchy +=  fac_gp * cauchy_al_3 ;

      //   NPK2 +=  fac_gp * PK2_al_be ;
      //   MPK2 +=  fac_gp * zeta *PK2_al_be ;
      //   MPK2_2 +=  fac_gp * zeta *PK2_al_be ;
      //   QPK2 +=  fac_gp *PK2_al_3;
      // }

      // return std::tuple(NCauchy,MCauchy,MCauchy2,QCauchy);
    }

    auto stressResultants(const ResultRequirementsType &req, const Dune::FieldVector<double, 2> &local,
                          ResultTypeMap<double> &result)const
    {
      // const auto& oneDRule = Dune::QuadratureRules<double,1>::rule(Dune::GeometryTypes::line,4);
      // const auto &thickness = this->fESettings.template request<double>("thickness");
      // using namespace Dune::DerivativeDirections;
      // using namespace Dune;
      // const auto [ displacementFunction, directorFunction,
      //     directorReferenceFunction]
      //     = this->template createFunctions<double>(req.getFERequirements());

      // Eigen::Matrix2d NCauchy = Eigen::Matrix2d::Zero();
      // Eigen::Matrix2d MCauchy = Eigen::Matrix2d::Zero();
      // Eigen::Vector2d QCauchy = Eigen::Vector2d::Zero();

      // Eigen::Matrix2d NPK2 = Eigen::Matrix2d::Zero();
      // Eigen::Matrix2d MPK2 = Eigen::Matrix2d::Zero();
      // Eigen::Vector2d QPK2 = Eigen::Vector2d::Zero();

      // KinematicVariables<double> kin;
      // for (auto& gP : oneDRule)
      // {
      //   const auto& gppos =gP.position();
      //   const auto& gpweight =gP.weight();

      //   const double zeta  = (2*gppos[0]-1)*thickness/2.0;
      //   const Dune::FieldVector<double,2> gp2DPos= local;
      //   const Dune::FieldVector<double,3> gp3DPos= {local[0],local[1],gppos[0]};

      //   const double fac_gp = gpweight * thickness * 0.5*2;
      //   // the first two fixes the change of the integration mapping from 0..1 to -1..1,
      //   // and the h/2 factor is the factor for the correct thickness

      //   const auto [res,detF,E11,zetaS,jS,JS,gp2DPosS]   = calculateStresses(req.getFERequirements(),gp3DPos,true);

      //   const auto [ cauchy,PK2,PK1] = res;
      //   Eigen::Matrix2d PK2_al_be= PK2.template block<2,2>(0,0);
      //   Eigen::Matrix2d cauchy_al_be= cauchy.template block<2,2>(0,0);

      //   const Eigen::Vector2d cauchy_al_3              = cauchy.template block<2,1>(0,2);
      //   const Eigen::Vector2d PK2_al_3              = PK2.template block<2,1>(0,2);

      //   NCauchy +=  fac_gp *  cauchy_al_be ;
      //   MCauchy +=  fac_gp *  zeta *cauchy_al_be ;
      //   QCauchy +=  fac_gp * cauchy_al_3 ;

      //   NPK2 +=  fac_gp * PK2_al_be ;
      //   MPK2 +=  fac_gp * zeta *PK2_al_be ;
      //   QPK2 +=  fac_gp *PK2_al_3;
      // }
      // typename ResultTypeMap<double>::ResultArray resultVector;

      // if (req.isResultRequested(ResultType::membraneForces)) {
      //   resultVector.resize(2, 2);
      //   resultVector=NCauchy;
      //   result.insertOrAssignResult(ResultType::membraneForces, resultVector);
      // }else if (req.isResultRequested(ResultType::bendingMoments)) {
      //   resultVector.resize(2, 2);
      //   resultVector=MCauchy;
      //   result.insertOrAssignResult(ResultType::bendingMoments, resultVector);
      // }else if (req.isResultRequested(ResultType::shearForces)) {
      //   resultVector.resize(2, 1);
      //   resultVector=QCauchy;
      //   result.insertOrAssignResult(ResultType::shearForces, resultVector);
      // }else  if (req.isResultRequested(ResultType::membraneForcesPK2)) {
      //   resultVector.resize(2, 2);
      //   resultVector=NPK2;
      //   result.insertOrAssignResult(ResultType::membraneForcesPK2, resultVector);
      // }else if (req.isResultRequested(ResultType::bendingMomentsPK2)) {
      //   resultVector.resize(2, 2);
      //   resultVector=MPK2;
      //   result.insertOrAssignResult(ResultType::bendingMomentsPK2, resultVector);
      // }else if (req.isResultRequested(ResultType::shearForcesPK2)) {
      //   resultVector.resize(2, 1);
      //   resultVector=QPK2;
      //   result.insertOrAssignResult(ResultType::shearForcesPK2, resultVector);
      // }
    }

    auto calculateInternalEnergy(const FERequirementType &par, const Dune::FieldVector<double, 2> &local) const  {
      const auto& oneDRule = Dune::QuadratureRules<double,1>::rule(Dune::GeometryTypes::line,4);

      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto [ displacementFunction, directorFunction,
          directorReferenceFunction]
          = this->template createFunctions<double>(par);
      const auto &lambda = par.getParameter(FEParameter::loadfactor);
      Eigen::Vector3d energy;
      energy.setZero();
      this->membraneStrain.pre(geo,displacementFunction);

      const auto &thickness_ = this->settings_.thickness;
      KinematicVariables<double> kin{};
//      for (auto& gP : oneDRule)
//      {
//        const auto& gppos =gP.position();
//        const auto& gpweight =gP.weight();
//        const double zeta  = (2*gppos[0]-1)*thickness_/2.0;
        const Dune::FieldVector<double,2> gp2DPos= local;
//        const Dune::FieldVector<double,3> gp3DPos= {local[0],local[1],gppos[0]};

        kin.t           = directorFunction.evaluate(gp2DPos,Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.t0          = directorReferenceFunction.evaluate(gp2DPos,Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.ud1andud2   = displacementFunction.evaluateDerivative(gp2DPos, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.A1andA2     = Dune::toEigen(geo.jacobianTransposed(gp2DPos)).transpose();
        kin.a1anda2     = kin.A1andA2+ kin.ud1andud2;
        kin.t0d1Andt0d2 = directorReferenceFunction.evaluateDerivative(gp2DPos, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));
        kin.td1Andtd2   = directorFunction.evaluateDerivative(gp2DPos, Dune::wrt(spatialAll),Dune::on(Dune::DerivativeDirections::referenceElement));

        auto [_,epsV,kappaV,gammaV]                = this->computeMaterialAndStrains(gp2DPos,0,geo,displacementFunction,directorFunction,directorReferenceFunction,kin);
        Eigen::Vector3<double> rhoV;
        rhoV<< 0.5*(kin.td1().squaredNorm()-kin.t0d1().squaredNorm()),0.5*(kin.td2().squaredNorm()-kin.t0d2().squaredNorm()),kin.td1().dot(kin.td2())-kin.t0d1().dot(kin.t0d2());
//        const auto strainsV= (epsV+ zeta*kappaV+secondOrderBending*zeta*zeta*rhoV).eval();

        Eigen::Matrix3d J;
//        const auto G1AndG2 = (kin.A1andA2+ zeta* kin.t0d1Andt0d2).eval();
        J<<kin.A1andA2.col(0),kin.A1andA2.col(1),kin.t0;

        const  Eigen::Matrix3d Jloc= orthonormalizeMatrixColumns(J);
        const Eigen::Matrix3d T= Jloc*J.inverse();

        Eigen::Matrix3d EV;
        EV<< epsV(0,0),epsV(2)/2, 0,
            epsV(2)/2, epsV(1), 0,
            0,0,0;
        EV= T.transpose()*EV*T;

        Eigen::Matrix3d KV;
        KV<< kappaV(0),(kappaV(2))/2, 0,
            (kappaV(2))/2, kappaV(1), 0,
            0,0,0;
        KV= T.transpose()*KV*T;

        Eigen::Matrix3d KV2;
        KV2<< rhoV(0),(rhoV(2))/2, 0,
            (rhoV(2))/2, rhoV(1), 0,
            0,0,0;
        KV2= T.transpose()*KV2*T;

        Eigen::Matrix3d GV;
        GV<< 0,0, gammaV(0)/2,
            0, 0, gammaV(1)/2,
            gammaV(0)/2,gammaV(1)/2,0;
//        GV= T.transpose()*GV*T;
        Eigen::Vector2d gaV;
        gaV<<GV(0,2)*2,GV(1,2)*2;
        Eigen::Matrix3d QCauchyBig;
        QCauchyBig.setZero();

//        const auto [res,detF,E11,zetaS,jS,JS,gp2DPosS]   = calculateStresses(par,gp3DPos,true);



//        const double fac_gp = gpweight * thickness_ * 0.5*2;



        const auto [NCauchy,MCauchy,MCauchy2,QCauchy] =calcStressResultants(par,local);
        QCauchyBig.template block<2,1>(0,2)=QCauchy;
        QCauchyBig.template block<1,2>(2,0)=QCauchy;
//        std::cout<<"QCauchyBig"<<std::endl;
//        std::cout<<QCauchyBig<<std::endl;
//        std::cout<<"GV"<<std::endl;
//        std::cout<<GV<<std::endl;
//        std::cout<<"gammaV"<<std::endl;
//        std::cout<<gammaV<<std::endl;
        const double membrane = 0.5*(NCauchy*EV.template block<2,2>(0,0)).trace();
        const double bendingE =0.5*(MCauchy*KV.template block<2,2>(0,0)+MCauchy2*KV2.template block<2,2>(0,0)).trace();
        const double shearE =0.5*(QCauchyBig*GV).trace();
        Eigen::Vector3d en;
        en<<membrane,bendingE,shearE;
        return en;
//        energy += (en)*fac_gp;
//      }

      return energy;
    }
  };




}  // namespace Ikarus