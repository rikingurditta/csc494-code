#include <igl/opengl/glfw/Viewer.h>

#include "get_decimated_sequence.h"
#include "flow.h"
#include "igl/readCSV.h"
#include "reinflate.h"
#include "closest_points.h"
#include "BAPN.h"
#include "helpers.h"

Eigen::MatrixXd make_3d_mat_from_2d(const Eigen::MatrixXd &V2d) {
    Eigen::MatrixXd V3 = Eigen::MatrixXd::Zero(V2d.rows(), 3);
    V3.block(0, 0, V2d.rows(), 2) = V2d;
    return V3;
}

void make_3d_mesh_from_2d(const Eigen::MatrixXd &V2d, const Eigen::MatrixXi &E2d,
                          Eigen::MatrixXd &V3, Eigen::MatrixXi &F3d) {
    // all z coordinates are 0
    V3 = make_3d_mat_from_2d(V2d);
    // each 2d mesh edge becomes a degenerate triangle
    F3d = Eigen::MatrixXi::Zero(E2d.rows(), 3);
    F3d.block(0, 0, E2d.rows(), 2) = E2d;
    F3d.col(2) = E2d.col(1);
}

void plot_grads(const Eigen::MatrixXd &G, const Eigen::MatrixXd &Vc, const Eigen::MatrixXd &Vf, double h,
                igl::opengl::glfw::Viewer &viewer, int index) {
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero((Vc.rows() + Vf.rows()) * 2, 3);
    Eigen::MatrixXi E(Vc.rows() + Vf.rows(), 2);
    Eigen::Vector3d red = Eigen::Vector3d(1, 0, 0);
    Eigen::Vector3d green = Eigen::Vector3d(0, 1, 0);
    Eigen::MatrixXd C(Vc.rows() + Vf.rows(), 3);
    C.block(0, 0, Vc.rows(), 3).rowwise() = red.transpose();
    C.block(Vc.rows(), 0, Vf.rows(), 3).rowwise() = green.transpose();
    V.block(0, 0, (Vc.rows() + Vf.rows()) * 2, 2) << Vc, Vf, Vc, Vf;
    V.block(Vc.rows() + Vf.rows(), 0, Vc.rows() + Vf.rows(), 2) += G * h;
    for (int i = 0; i < Vc.rows() + Vf.rows(); i++) {
        E(i, 0) = i;
        E(i, 1) = i + Vc.rows() + Vf.rows();
    }
    viewer.data(index).set_edges(V, E, C);
}

void nested_cages() {
    Eigen::MatrixXd Vf, Vc;
    Eigen::MatrixXi Ef, Ec;
    igl::readCSV("../data/octagon_hole_hand_V.csv", Vf);
    igl::readCSV("../data/octagon_hole_hand_E.csv", Ef);
    igl::readCSV("../data/square_hole_hand_V.csv", Vc);
    igl::readCSV("../data/square_hole_hand_E.csv", Ec);


    int max_flow_meshes = 24;
    std::vector<Eigen::MatrixXd> Vc_flow;
    std::vector<Eigen::MatrixXd> Vf_flow;
    std::vector<Eigen::MatrixXd> reinflation_grads;
    int num_flow_meshes = flow(Vc, Ec, Vf, Ef, max_flow_meshes, Vf_flow, 0.01);
    for (int i = 0; i < num_flow_meshes + 1; i++) {
        Vc_flow.emplace_back(Vc);
        reinflation_grads.emplace_back(Eigen::MatrixXd::Zero(Vc.rows() + Vf.rows(), 2));
    }
    reinflation_grads.pop_back();
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
    viewer.append_mesh();
    viewer.append_mesh();
    viewer.append_mesh();
    viewer.append_mesh();
    viewer.data(0).set_mesh(Vf_viewer, Ff_viewer);
    viewer.data(1).set_mesh(Vc_viewer, Fc_viewer);
    viewer.data(4).set_mesh(Vf_viewer, Ff_viewer);
    viewer.data(4).set_colors(Eigen::MatrixXd::Ones(Vf_viewer.rows(), 3));

    Eigen::MatrixXd C(Vf.rows(), 2);
    query_closest_points(Vf, Vc, Ec, C);
    Eigen::MatrixXd V_closest_points_viewer = make_3d_mat_from_2d(C);
    Eigen::MatrixXd closest_points_colours = Eigen::MatrixXd::Ones(C.rows(), 3);
    viewer.data(2).set_points(V_closest_points_viewer, closest_points_colours);
    viewer.data(2).point_size = 5;

    Eigen::MatrixXd outer_square(Vc.rows() / 2, 2);
    Eigen::MatrixXd V_outer_square_viewer = Eigen::MatrixXd::Zero(Vc.rows() / 2, 3);
    V_outer_square_viewer.block(0, 0, V_outer_square_viewer.rows(), 2)
            = Vc.block(0, 0, Vc.rows() / 2, 2);
    Eigen::MatrixXd outer_square_colours = Eigen::MatrixXd::Zero(V_outer_square_viewer.rows(), 3);
    outer_square_colours.block(0, 0, outer_square_colours.rows(), 1) = Eigen::MatrixXd::Ones(
            outer_square_colours.rows(), 1);
    viewer.data(3).set_points(V_outer_square_viewer, outer_square_colours);
    viewer.data(3).point_size = 5;

    Eigen::VectorXd x = Eigen::VectorXd::Zero(Vc.rows() * 2 + Vf.rows() * 2);
    x.segment(0, Vc.rows() * 2) = flatten(Vc);
    x.segment(Vc.rows() * 2, Vf.rows() * 2) = flatten(Vf_flow[num_flow_meshes]);
    Eigen::VectorXd dir;
    for (int i = num_flow_meshes; i > 0; i--) {
        std::cout << "i: " << i << "\n";
        Eigen::MatrixXd Vf_next = Vf_flow[i - 1];
        Eigen::MatrixXd Vf_curr = unflatten(x, 2).block(Vc.rows(), 0, Vf.rows(), 2);
        // descent iterations
        for (int j = 0; j < 100; j++) {
            auto [E, g, H] = total_energy(Vc, Ec, Vf_curr, Vf_next, Ef, x, 0.1);

            dir = -g;
            // could check norm of grad to see if converged/can't step
            double alpha = 0.1;
            double beta = 0.5;
            double t = 0.1;
            Eigen::VectorXd x_temp;
            double E_temp;
            do {
                t *= beta;
                x_temp = x + t * dir;
                auto [E_, g_, H_] = total_energy(Vc, Ec, Vf_curr, Vf_next, Ef, x_temp, 0.1);
                E_temp = E_;
            } while (E_temp > E + alpha * t * g.dot(dir) and t > 1e-15);
            if (t <= 1e-15) {
                // don't need to step if converged
                // need to throw error if not converged - can't step
                // TODO
            }
            x = x_temp;
        }
        reinflation_grads.emplace_back(unflatten(dir, 2));
        Eigen::MatrixXd meshes = unflatten(x, 2);
        Vc_flow.emplace_back(meshes.block(0, 0, Vc.rows(), 2));
        Vf_flow.emplace_back(meshes.block(Vc.rows(), 0, Vf.rows(), 2));
    }
    reinflation_grads.emplace_back(Eigen::MatrixXd::Zero(Vc.rows() + Vf.rows(), 2));


    std::cout << "press a/s to reverse/progress flow\npress d to go to end of flow/beginning of reinflation\n";
    long m = 0;
    viewer.callback_key_pressed =
            [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod) {
                switch (key) {
                    case 'a':
                        m--;
                        break;
                    case 's':
                        m++;
                        break;
                    case 'd':
                        m = num_flow_meshes;
                        break;
                    default:
                        return false;
                }
                if (m < 0)
                    m = 0;
                else if (m >= Vf_flow.size())
                    m = Vf_flow.size() - 1;

                std::cout << "m: " << m << "\n";
                Vf_viewer.block(0, 0, Vf.rows(), 2) = Vf_flow[m];
                viewer.data(0).set_mesh(Vf_viewer, Ff_viewer);
                Vc_viewer.block(0, 0, Vc.rows(), 2) = Vc_flow[m];
                viewer.data(1).set_mesh(Vc_viewer, Fc_viewer);
                Eigen::MatrixXd Vf_target_viewer = Vf_viewer;
                if (m >= num_flow_meshes) {
                    Vf_target_viewer.block(0, 0, Vf.rows(), 2)
                            = Vf_flow[num_flow_meshes * 2 - m];
                }
                viewer.data(4).set_mesh(Vf_target_viewer, Ff_viewer);

                query_closest_points(Vf_flow[m], Vc_flow[m], Ec, C);
                V_closest_points_viewer.block(0, 0, V_closest_points_viewer.rows(), 2) = C;
                viewer.data(2).set_points(V_closest_points_viewer, closest_points_colours);

                V_outer_square_viewer.block(0, 0, V_outer_square_viewer.rows(), 2) = Vc_viewer.block(0, 0,
                                                                                                     V_outer_square_viewer.rows(),
                                                                                                     2);
                viewer.data(3).set_points(V_outer_square_viewer, outer_square_colours);

                plot_grads(reinflation_grads[m], Vc_flow[m], Vf_flow[m], 100, viewer, 5);
                return true;
            };
    viewer.launch();
}

int main(int argc, char *argv[]) {
    nested_cages();
}