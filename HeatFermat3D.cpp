#include <array>
#include "Spirals.h"


#include "SpiralDt.h"


int main(int argc, char *argv[]) {
    TriMesh  m;
    IterLine iter_line;
    OpenMesh::IO::read_mesh(m, argv[1]);


    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    std::tie(V,F) = fv_indices(m);
    Eigen::VectorXd d = distance_field_from_heat(V, F, 0);

    isoline_cut(m, iter_line, d, 200);

    std::ofstream obj("../../isoline.obj");
    obj << iter_line;
    obj.close();

    obj.open("../../spiral.obj");
    auto regions = iter_line.region_extract();
    auto [fs_start, fs_end] = iter_line.make_fermat_spiral(regions);
    iter_line.print_fermat_spiral(obj, fs_start);
    obj.close();

    m.request_face_normals();
    m.request_vertex_normals();
    m.update_face_normals();
    m.update_vertex_normals();
    for(auto vit = m.vertices_begin(); vit != m.vertices_end(); vit++) {
        m.point(*vit) -= 0.005*m.normal(*vit);
    }
    OpenMesh::IO::write_mesh(m, "../../out.off");
}
