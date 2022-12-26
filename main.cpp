#include <igl/opengl/glfw/Viewer.h>
#include <chrono>

using Clock = std::chrono::system_clock;


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

void plot_grads(const Eigen::MatrixXd &G, const Eigen::MatrixXd &Vc, const Eigen::MatrixXd &Vf,
                igl::opengl::glfw::Viewer &viewer, int index) {
    return;
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero((Vc.rows() + Vf.rows()) * 2, 3);
    Eigen::MatrixXi E(Vc.rows() + Vf.rows(), 2);
    Eigen::RowVector3d red = Eigen::Vector3d(0.8, 0, 0);
    Eigen::RowVector3d green = Eigen::Vector3d(0, 0.8, 0);
    Eigen::MatrixXd C(Vc.rows() + Vf.rows(), 3);
    C.block(0, 0, Vc.rows(), 3).rowwise() = red;
    C.block(Vc.rows(), 0, Vf.rows(), 3).rowwise() = green;
    V.block(0, 0, (Vc.rows() + Vf.rows()) * 2, 2) << Vc, Vf, Vc, Vf;
    V.block(Vc.rows() + Vf.rows(), 0, Vc.rows() + Vf.rows(), 2) += G.rowwise().normalized() * 0.25;
    for (int i = 0; i < Vc.rows() + Vf.rows(); i++) {
        E(i, 0) = i;
        E(i, 1) = i + Vc.rows() + Vf.rows();
    }
    viewer.data(index).set_edges(V, E, C);
}

void nested_cages(int argc, char *argv[]) {
    if (argc < 3) {
        std::cout << "cages <fine mesh basename> <coarse mesh basename>\n";
        return;
    }
    Eigen::MatrixXd Vf, Vc;
    Eigen::MatrixXi Ef, Ec;
    igl::readCSV("../data/" + std::string(argv[1]) + "_V.csv", Vf);
    igl::readCSV("../data/" + std::string(argv[1]) + "_E.csv", Ef);
    igl::readCSV("../data/" + std::string(argv[2]) + "_V.csv", Vc);
    igl::readCSV("../data/" + std::string(argv[2]) + "_E.csv", Ec);

    auto start = Clock::now();

    int max_flow_meshes = 100;
    std::vector<Eigen::MatrixXd> Vc_flow;
    std::vector<Eigen::MatrixXd> Vf_flow;
    std::vector<Eigen::MatrixXd> grads;

    // flow step
    start = Clock::now();
    int num_flow_meshes = flow(Vc, Ec, Vf, Ef, max_flow_meshes, Vf_flow, grads, 0.01);
    auto flow_duration = std::chrono::duration<double>(Clock::now() - start);

    for (int i = 0; i < num_flow_meshes + 1; i++) {
        Vc_flow.emplace_back(Vc);
    }

    // reinflate step
    start = Clock::now();
    Eigen::VectorXd x = Eigen::VectorXd::Zero(Vc.rows() * 2 + Vf.rows() * 2);
    x.segment(0, Vc.rows() * 2) = flatten(Vc);
    x.segment(Vc.rows() * 2, Vf.rows() * 2) = flatten(Vf_flow[num_flow_meshes]);
    Eigen::VectorXd dir;

    for (int i = num_flow_meshes; i > 0; i--) {
//        std::cout << "i: " << i << "\n";
        Eigen::MatrixXd Vf_next = Vf_flow[i - 1];
        Eigen::MatrixXd Vf_curr = unflatten(x, 2).block(Vc.rows(), 0, Vf.rows(), 2);
        // descent iterations
        for (int j = 0; j < 100; j++) {
//            std::cout << "j: " << j << "\n";
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
        grads.emplace_back(unflatten(dir, 2));
        Eigen::MatrixXd meshes = unflatten(x, 2);
        Vc_flow.emplace_back(meshes.block(0, 0, Vc.rows(), 2));
        Vf_flow.emplace_back(meshes.block(Vc.rows(), 0, Vf.rows(), 2));
    }
    auto reinflate_duration = std::chrono::duration<double>(Clock::now() - start);

    grads.emplace_back(Eigen::MatrixXd::Zero(Vc.rows() + Vf.rows(), 2));


    // make 3d versions of 2d mesh to display in viewer
    Eigen::MatrixXd Vf_viewer, Vc_viewer;
    Eigen::MatrixXi Ff_viewer, Fc_viewer;
    make_3d_mesh_from_2d(Vf, Ef, Vf_viewer, Ff_viewer);
    make_3d_mesh_from_2d(Vc, Ec, Vc_viewer, Fc_viewer);
    Eigen::MatrixXd Ef_colours(Ef.rows(), 3);
    Ef_colours.rowwise() = Eigen::RowVector3d(0, 0.5, 0);
    Eigen::MatrixXd Ec_colours(Ec.rows(), 3);
    Ec_colours.rowwise() = Eigen::RowVector3d(0.5, 0, 0);

    igl::opengl::glfw::Viewer viewer;
    viewer.core().background_color.setOnes();
    viewer.append_mesh();
    viewer.append_mesh();
    viewer.append_mesh();
    viewer.append_mesh();
    viewer.data(0).set_edges(Vf_viewer, Ef, Ef_colours);
    viewer.data(1).set_edges(Vc_viewer, Ec, Ec_colours);
    viewer.data(3).set_mesh(Vf_viewer, Ff_viewer);
    viewer.data(3).set_colors(Eigen::MatrixXd::Ones(Vf_viewer.rows(), 3));

    Eigen::MatrixXd C(Vf.rows(), 2);
    query_closest_points(Vf, Vc, Ec, C);
    Eigen::MatrixXd V_closest_points_viewer = make_3d_mat_from_2d(C);
    Eigen::MatrixXd closest_points_colours = Eigen::MatrixXd::Ones(C.rows(), 3);
    viewer.data(2).set_points(V_closest_points_viewer, closest_points_colours);
    viewer.data(2).point_size = 5;
    plot_grads(grads[0], Vc_flow[0], Vf_flow[0], viewer, 5);

    std::cout << "flow time: " << flow_duration.count() << " s\t"
              << "reinflate time: " << reinflate_duration.count() << " s\n";

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
                viewer.data(0).set_edges(Vf_viewer, Ef, Ef_colours);
                Vc_viewer.block(0, 0, Vc.rows(), 2) = Vc_flow[m];
                viewer.data(1).set_edges(Vc_viewer, Ec, Ec_colours);
                Eigen::MatrixXd Vf_target_viewer = Vf_viewer;
                if (m >= num_flow_meshes) {
                    Vf_target_viewer.block(0, 0, Vf.rows(), 2)
                            = Vf_flow[num_flow_meshes * 2 - m];
                }
                viewer.data(3).set_mesh(Vf_target_viewer, Ff_viewer);

                query_closest_points(Vf_flow[m], Vc_flow[m], Ec, C);
                V_closest_points_viewer.block(0, 0, V_closest_points_viewer.rows(), 2) = C;
                viewer.data(2).set_points(V_closest_points_viewer, closest_points_colours);

                plot_grads(grads[m], Vc_flow[m], Vf_flow[m], viewer, 4);
                return true;
            };
    viewer.launch();
}

int main(int argc, char *argv[]) {
    nested_cages(argc, argv);
}