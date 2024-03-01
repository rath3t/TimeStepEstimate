// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <dune/common/float_cmp.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/eigen.h>
#include <dune/python/pybind11/operators.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/alugrid/grid.hh>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/python/finiteelements/materials/material.hh>
#include <ikarus/utils/basis.hh>

#include "rmstiffness.hh"
#include "resultantrmshell.hh"



int add(int i, int j) {
    return i + j;
}

class Matrix {
public:
    Matrix(size_t rows, size_t cols) : m_rows(rows), m_cols(cols) {
        m_data = new float[rows*cols];
    }
    float *data() { return m_data; }
    size_t rows() const { return m_rows; }
    size_t cols() const { return m_cols; }

    float &operator()(const std::array<int,2>& ar) { return m_data[ar[0]*m_cols+ar[1]]; }
private:
    size_t m_rows, m_cols;
    float *m_data;
};

PYBIND11_MODULE(_timestepestimate, m) {
      m.doc() = "Python bindings for time step estimates"; // optional module docstring

  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace Ikarus;

  py::class_<Matrix>(m, "Matrix", py::buffer_protocol())
  .def(py::init<int,int>())
  .def("__getitem__", &Matrix::operator())
  .def("__setitem__", [](Matrix& self, const std::array<int,2>& ar, float val) { self(ar) = val; })
   .def_buffer([](Matrix &m) -> py::buffer_info {
        return py::buffer_info(
            m.data(),                               /* Pointer to buffer */
            sizeof(float),                          /* Size of one scalar */
            py::format_descriptor<float>::format(), /* Python struct-style format descriptor */
            2,                                      /* Number of dimensions */
            { m.rows(), m.cols() },                 /* Buffer dimensions */
            { sizeof(float) * m.cols(),             /* Strides (in bytes) for each index */
              sizeof(float) }
        );
    });


   m.def("computeStiffnessMatrix", [](py::dict settings,const Eigen::Ref<Eigen::Matrix3Xd>& referenceMidSurface,const Eigen::Ref<Eigen::Matrix3Xd>& referenceDirectors, const Eigen::Ref<Eigen::Matrix3Xd>& displacements, Eigen::Ref<Eigen::Matrix3Xd>& directors) {
std::cout<<referenceMidSurface<<std::endl;
std::cout<<referenceDirectors<<std::endl;
std::cout<<displacements<<std::endl;
std::cout<<directors<<std::endl;
         using Grid = Dune::ALUGrid<2, 3, Dune::cube, Dune::nonconforming>;
        using namespace Dune::Functions::BasisFactory;

    using GridView = Grid::LeafGridView;
    using Basis =  decltype(Ikarus::makeBasis(
      std::declval<GridView>(), composite(power<3>(lagrange<1>(), BlockedInterleaved()),
                          power<2>(lagrange<1>(), BlockedInterleaved()), BlockedLexicographic{})));
    
      using MidSurfaceVector   = Dune::BlockVector<Dune::RealTuple<double, 3>>;
  using DirectorVector     = Dune::BlockVector<Dune::UnitVector<double, 3>>;
  using MultiTypeVector    = Dune::TupleVector<MidSurfaceVector, DirectorVector>;

  Dune::GridFactory<Grid> gridFactory;
  for (size_t i = 0; i < referenceMidSurface.cols(); i++)
  {
   gridFactory.insertVertex({referenceMidSurface(0, i), referenceMidSurface(1, i), referenceMidSurface(2, i)});
  }

  
  
  gridFactory.insertElement(Dune::GeometryTypes::cube(2),{0,1,2,3});
  auto grid = gridFactory.createGrid();
  auto gridView = grid->leafGridView();
  using namespace Dune::Functions::BasisFactory;

auto scalarMidSurfBasis = lagrange<1>();
auto scalaarDirectorBasis = lagrange<1>();
   auto basis = Ikarus::makeBasis(
      gridView, composite(power<3>(scalarMidSurfBasis, BlockedInterleaved()),
                          power<2>(scalaarDirectorBasis, BlockedInterleaved()), BlockedLexicographic{}));

    MidSurfaceVector displacementsBlockedRef(basis.untouched().size({Dune::Indices::_0}));

    MidSurfaceVector displacementsBlocked(basis.untouched().size({Dune::Indices::_0}));


  for (int i =0;auto& dispSingle : displacementsBlocked) {
  dispSingle.setValue(displacements.col(i++).eval());
  }

    for (auto& dispSingle : displacementsBlockedRef) {
  dispSingle.setValue(Eigen::Vector<double, 3>::Zero());
  }


  DirectorVector dBlocked(basis.untouched().size({Dune::Indices::_1}));
  for (int i =0;auto& dsingle : dBlocked) {
  dsingle.setValue(directors.col(i++));
  }

    DirectorVector dBlockedRef(basis.untouched().size({Dune::Indices::_1}));
  for (int i =0;auto& dsingle : dBlockedRef) {
  dsingle.setValue(referenceDirectors.col(i++));
  }

  const MultiTypeVector x0(displacementsBlockedRef, dBlockedRef);
   MultiTypeVector x(displacementsBlocked, dBlocked);

using FE = Ikarus::StressBasedShellRM<Ikarus::NonLinearRMshell<Basis,true>>;

RMSettings rmSettings;
rmSettings.thickness = settings["thickness"].cast<double>();
rmSettings.Emodul = settings["youngsModulus"].cast<double>();
rmSettings.nu = settings["nu"].cast<double>();
std::cout<<"x0"<<x0<<std::endl;
std::cout<<"x"<<x<<std::endl;

double load =0.0;
auto fe = FE(basis,*gridView.begin<0>(),x0,rmSettings);
  auto req = Ikarus::FERequirements<std::reference_wrapper<MultiTypeVector>>().addAffordance(
      Ikarus::AffordanceCollections::elastoStatics);
        req.insertGlobalSolution(Ikarus::FESolutions::noSolution, x)
        .insertParameter(Ikarus::FEParameter::loadfactor,load);

        Eigen::MatrixXd stiffnessMatrix;
        stiffnessMatrix.setZero(fe.size(),fe.size());
fe.calculateMatrix(req,stiffnessMatrix);
std::cout<<stiffnessMatrix<<std::endl;
     return stiffnessMatrix;
     });




    m.def("add", &add, "A function which adds two numbers",
      py::arg("i"), py::arg("j"));

}
