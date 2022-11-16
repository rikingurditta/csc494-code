#include <igl/opengl/glfw/Viewer.h>

#include "get_decimated_sequence.h"
#include "flow.h"
#include "igl/readCSV.h"
#include "reinflate.h"
#include "closest_points.h"
#include "BAPN.h"
#include "helpers.h"

void make_3d_mesh_from_2d(const Eigen::MatrixXd &V2d, const Eigen::MatrixXi &E2d,
                          Eigen::MatrixXd &V3, Eigen::MatrixXi &F3d)
{
    // all z coordinates are 0
    V3 = Eigen::MatrixXd::Zero(V2d.rows(), 3);
    V3.block(0, 0, V2d.rows(), 2) = V2d;
    // each 2d mesh edge becomes a degenerate triangle
    F3d = Eigen::MatrixXi::Zero(E2d.rows(), 3);
    F3d.block(0, 0, E2d.rows(), 2) = E2d;
    F3d.col(2) = E2d.col(1);
}

void nested_cages()
{
    Eigen::MatrixXd Vf, Vc;
    Eigen::MatrixXi Ef, Ec;
    igl::readCSV("../data/2d_idk_V.csv", Vf);
    igl::readCSV("../data/2d_idk_E.csv", Ef);
    igl::readCSV("../data/2d_idk_dec1_V.csv", Vc);
    igl::readCSV("../data/2d_idk_dec1_E.csv", Ec);

    int max_flow_meshes = 24;
    std::vector<Eigen::MatrixXd> Vf_flow;
    double h = 0.05;
    int num_flow_meshes = flow(Vc, Ec, Vf, Ef, max_flow_meshes, Vf_flow, h);
    Eigen::MatrixXd Vc_inflated;
//    reinflate(Vc, Ec, Vf, Ef, Vf_flow, Vc_inflated);

    // make 3d versions of 2d mesh to display in viewer
    Eigen::MatrixXd Vf_viewer, Vc_viewer;
    Eigen::MatrixXi Ff_viewer, Fc_viewer;
    make_3d_mesh_from_2d(Vf, Ef, Vf_viewer, Ff_viewer);
    make_3d_mesh_from_2d(Vc, Ec, Vc_viewer, Fc_viewer);

    igl::opengl::glfw::Viewer viewer;
    viewer.append_mesh();
    viewer.append_mesh();
    viewer.data(0).set_mesh(Vf_viewer, Ff_viewer);
    viewer.data(1).set_mesh(Vc_viewer, Fc_viewer);

    Eigen::MatrixXd C(Vf.rows(), 2);
    query_closest_points(Vf, Vc, Ec, C);
    Eigen::MatrixXd V_closest_points_viewer = Eigen::MatrixXd::Zero(C.rows(), 3);
    V_closest_points_viewer.block(0, 0, V_closest_points_viewer.rows(), 2) = C;
    Eigen::MatrixXd closest_points_colours = Eigen::MatrixXd::Ones(C.rows(), 3);
    viewer.data(2).set_points(V_closest_points_viewer, closest_points_colours);
    viewer.data(2).point_size = 5;

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
                else if (m > num_flow_meshes)
                    m = num_flow_meshes;

                std::cout << m << "\n";
                Vf_viewer.block(0, 0, Vf.rows(), 2) = Vf_flow[m];
                viewer.data(0).set_mesh(Vf_viewer, Ff_viewer);
                query_closest_points(Vf_flow[m], Vc, Ec, C);
                V_closest_points_viewer.block(0, 0, V_closest_points_viewer.rows(), 2) = C;
                viewer.data(2).set_points(V_closest_points_viewer, closest_points_colours);
                return true;
            };
    viewer.launch();
}

void test_bapn()
{
    Eigen::MatrixXd Vf, Vc;
    Eigen::MatrixXi Ef, Ec;
    igl::readCSV("../data/square_hole_hand_V.csv", Vf);
    igl::readCSV("../data/square_hole_hand_E.csv", Ef);
    igl::readCSV("../data/2d_idk_dec1_V.csv", Vc);
    igl::readCSV("../data/2d_idk_dec1_E.csv", Ec);

    Eigen::VectorXd x_query = flatten(Vf);
    std::cout << Vf << "\n";
    std::cout << Ef << "\n";

//    double E = nested_cages_energy(Vf, Ef, x_query);
    auto E = nested_cages_energy(Vf, Ef, x_query);
    std::cout << "Energy: " << std::get<0>(E) << "\n";
}

int main(int argc, char *argv[])
{
//    nested_cages();
    test_bapn();
}