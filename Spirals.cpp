//
// Created by hotpot on 2020/3/24.
//

#include "Spirals.h"

#include <igl/per_vertex_normals.h>
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <igl/writeDMAT.h>
#include <igl/opengl/glfw/Viewer.h>

#include <iostream>


std::tuple<Eigen::MatrixXd, Eigen::MatrixXi> fv_indices(TriMesh& m) {
    Eigen::MatrixXd V(m.n_vertices(), 3);
    Eigen::MatrixXi F(m.n_faces(), 3);
    for(int i = 0; i < V.rows(); i++)
        V.row(i) = *(m.points()+i);
    for(int i = 0; i < F.rows(); i++) {
        auto fv_it = m.fv_cwbegin(m.face_handle(i));
        F.row(i) << (fv_it++)->idx(), (fv_it++)->idx(), (fv_it++)->idx();
    }
    return {V, F};
}

Eigen::VectorXd distance_field_from_heat(Eigen::MatrixXd& V, Eigen::MatrixXi& F, int heat_source_idx) {
    igl::HeatGeodesicsData<double> data;
    if(!igl::heat_geodesics_precompute(V,F,0.1,data)) {
        std::cerr<<"Error: heat_geodesics_precompute failed."<<std::endl;
        exit(EXIT_FAILURE);
    };

    Eigen::VectorXd D = Eigen::VectorXd::Zero(data.Grad.cols());
    igl::heat_geodesics_solve(data,(Eigen::VectorXi(1,1)<<heat_source_idx).finished(),D);
    return D;
}



void isoline_cut(TriMesh& m, IterLine& line, Eigen::VectorXd& D, int n) {
    const double min = D.minCoeff(), max = D.maxCoeff();
    if (n == -1) {
        double avg_d_dis = std::accumulate(m.halfedges_begin(), m.halfedges_end(), 0.0f,
           [&](const double acc, const TriMesh::HalfedgeHandle eh){
               return acc + abs(D((m.from_vertex_handle(eh).idx())-D((m.to_vertex_handle(eh).idx()))));
           })/m.n_halfedges();
        n = (max-min)/avg_d_dis;
    }
    n += 1;

    std::vector<bool> visited(m.n_edges());
    for(int iso_i = 1; iso_i < n; iso_i++) {
        double iso_v = min + iso_i * (max-min)/n;
        std::fill(visited.begin(), visited.end(), false);
        for(auto eit = m.edges_begin(); eit != m.edges_end(); eit++) {
            auto _heh = m.halfedge_handle(*eit, 0);
            if (visited[eit->idx()] || D(m.from_vertex_handle(_heh).idx()) > iso_v
                                     || D(m.to_vertex_handle(_heh).idx()) < iso_v) continue;
            std::vector<Eigen::RowVector3d> isoline;
            TriMesh::HalfedgeHandle cur_heh = _heh;
            do {
                visited[m.edge_handle(cur_heh).idx()] = true;
                auto _fh = m.from_vertex_handle(cur_heh), _th = m.to_vertex_handle(cur_heh);
                double alpha = (iso_v - D(_fh.idx())) / (D(_th.idx()) - D(_fh.idx()));
                isoline.push_back(m.point(_fh) * (1 - alpha) + m.point(_th) * alpha);

                auto _oh = m.opposite_vh(cur_heh);
                double alpha_next = (iso_v - D(_th.idx())) / (D(_oh.idx()) - D(_th.idx()));
                double alpha_prev = (iso_v - D(_oh.idx())) / (D(_fh.idx()) - D(_oh.idx()));
                if (alpha_next >= 0 and alpha_next <= 1) {
                    cur_heh = m.next_halfedge_handle(cur_heh);
                } else if (alpha_prev >= 0 and alpha_prev <= 1) {
                    cur_heh = m.prev_halfedge_handle(cur_heh);
                } else {
                    std::cerr << "unexpected value catched " << alpha_next << ' ' << alpha_prev << std::endl;
                }
                cur_heh = m.opposite_halfedge_handle(cur_heh);
            } while (cur_heh != _heh);
            line.add_hierarchy(line.add_isoline(isoline), iso_i);
        }
    }
}

