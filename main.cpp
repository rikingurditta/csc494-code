#include <igl/opengl/glfw/Viewer.h>

#include "get_decimated_sequence.h"
#include "flow.h"
#include "reinflate.h"

int main(int argc, char *argv[])
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOBJ("../data/stanford-bunny.obj", V, F);
    std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> decimated_meshes;

    if (not get_decimated_sequence(V, F, 2, F.rows() * 0.01, decimated_meshes))
        return 1;

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(decimated_meshes[0].first, decimated_meshes[0].second);
    viewer.data().set_face_based(true);
    viewer.launch();
}
