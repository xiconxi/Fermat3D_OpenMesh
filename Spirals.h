//
// Created by hotpot on 2020/3/24.
//

#ifndef VISMAGNUM_SPIRALS_H
#define VISMAGNUM_SPIRALS_H

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/EigenVectorT.hh>
#include <igl/heat_geodesics.h>

#include "SpiralDt.h"

struct EigenTraits : OpenMesh::DefaultTraits {
    using Point = Eigen::Vector3d;
    using Normal = Eigen::Vector3d;
};
using TriMesh = OpenMesh::TriMesh_ArrayKernelT<EigenTraits>;

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi> fv_indices(TriMesh& m);

Eigen::VectorXd distance_field_from_heat(Eigen::MatrixXd& V, Eigen::MatrixXi& F, int heat_source_idx);

void isoline_cut(TriMesh& m, IterLine& line, Eigen::VectorXd& D, int n = -1);


#endif //VISMAGNUM_SPIRALS_H
