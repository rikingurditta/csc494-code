#include <igl/opengl/glfw/Viewer.h>

#include "get_decimated_sequence.h"
#include "flow.h"
#include "igl/readCSV.h"
#include "reinflate.h"

int main(int argc, char *argv[]) {
    Eigen::MatrixXd Vf, Vc;
    Eigen::MatrixXi Ef, Ec;
    igl::readCSV("../data/2d_torus_V.csv", Vf);
    igl::readCSV("../data/2d_torus_E.csv", Ef);
    igl::readCSV("../data/2d_torus_dec1_V.csv", Vc);
    igl::readCSV("../data/2d_torus_dec1_E.csv", Ec);

    std::vector<Eigen::MatrixXd> Vf_flow;
    flow(Vc, Ec, Vf, Ef, Vf_flow);
    Eigen::MatrixXd Vc_inflated;
//    reinflate(Vc, Ec, Vf, Ef, Vf_flow, Vc_inflated);

    // make 3d versions of 2d mesh to display in viewer
    // all z coordinates are 0
    Eigen:: MatrixXd V_viewer = Eigen::MatrixXd::Zero(Vf.rows(), 3);
    V_viewer.block(0, 0, Vf.rows(), 2) = Vf;
    // each 2d mesh edge becomes a degenerate triangle
    Eigen::MatrixXi F_viewer(Ef.rows(), 3);
    F_viewer.block(0, 0, F_viewer.rows(), 2) = Ef;
    F_viewer.col(2) = Ef.col(1);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V_viewer, F_viewer);

    std::cout << "press a/s to reverse/progress flow\n";
    int m = 0;
    viewer.callback_key_pressed =
            [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod) {
                switch (key) {
                    case 'a':
                        m--;
                        break;
                    case 's':
                        m++;
                        break;
                    default:
                        return false;
                }
                if (m < 0)
                    m = 0;
                else if (m > 25)
                    m = 25;

                std::cout << m << "\n";
                V_viewer.block(0, 0, Vf.rows(), 2) = Vf_flow[m];
                viewer.data().set_mesh(V_viewer, F_viewer);
                return true;
            };
    viewer.launch();
}
