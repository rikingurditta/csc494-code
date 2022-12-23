#include <iostream>
#include "flow.h"
#include "grad_sdf.h"
#include "igl/adjacency_list.h"
#include "closest_points.h"
#include "helpers.h"

void grad_unsigned_distance(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::RowVector2d &p,
                            Eigen::RowVector2d &g)
{
    Eigen::RowVector2d closest_point;
    closest_point_on_mesh(V, E, p, closest_point);
    g = (p - closest_point).normalized();
}

// doing quadrature by sampling points linearly along edges and weighting them by 1 / edge length and
// quadratic falloff from centre point
void quadrature_grad_dist(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E,
                          const Eigen::RowVector2d &p, const Eigen::RowVector2d &p_e0, const Eigen::RowVector2d &p_e1,
                          int num_points,
                          Eigen::RowVector2d &g)
{
    g = Eigen::RowVector2d::Zero();
    double l0 = (p - p_e0).norm();
    double l1 = (p - p_e1).norm();
    Eigen::RowVector2d g0;
    Eigen::RowVector2d g1;
//    grad_unsigned_distance(V, E, (p + p_e0) / 2, g0);
//    grad_unsigned_distance(V, E, (p + p_e1) / 2, g1);
//    g = (g0 / l0 + g1 / l1).normalized();
//    return;
    for (int i = 0; i <= num_points; i++) {
        // t is fraction of distance across edge
        double t = i * 1. / num_points;
        // get closest point direction sample at t% along either edge
        grad_unsigned_distance(V, E, (1. - t) * p + t * p_e0, g0);
        grad_unsigned_distance(V, E, (1. - t) * p + t * p_e1, g1);
        // add to total gradient with weight t^2 (also weight by edge length)
        g += pow(1. - t, 2) * (g0 / l0 + g1 / l1);
    }
    g.normalize();
}

bool
query_mesh_inside(const Eigen::MatrixXd &V_outside, const Eigen::MatrixXi &E_outside, const Eigen::MatrixXd &V_inside) {
    // TODO: this only checks if all vertices are inside, not if entire mesh is inside
    for (int v = 0; v < V_inside.rows(); v++) {
        if (!query_point_inside(V_outside, E_outside, V_inside.row(v))) {
//            std::cout << V_inside.row(v) << "\n";
            return false;
        }
    }
    return true;
}

int flow(const Eigen::MatrixXd &Vc, const Eigen::MatrixXi &Ec, const Eigen::MatrixXd &Vf, const Eigen::MatrixXi &Ef,
         int max_num_meshes, std::vector<Eigen::MatrixXd> &Vf_flow, double h) {
    int iteration = 0;
    Vf_flow.emplace_back(Vf);
    Eigen::MatrixXd V_curr = Vf;
    std::vector<std::vector<int>> A;
    igl::adjacency_list(Ef, A);
    while (!query_mesh_inside(Vc, Ec, V_curr) && iteration < max_num_meshes) {
        // g = ∇ sdf
        Eigen::MatrixXd g = Eigen::MatrixXd::Zero(Vf.rows(), Vf.cols());
        for (int v = 0; v < V_curr.rows(); v++) {
            Eigen::RowVector2d g_curr;
            quadrature_grad_dist(Vc, Ec,
                                 V_curr.row(v), V_curr.row(A[v][0]), V_curr.row(A[v][1]),
                                 4,
                                 g_curr);
//            grad_unsigned_distance(Vc, Ec, V_curr.row(v), g_curr);

            // calculate winding number
            int w = query_winding_number(Vc, Ec, V_curr.row(v));
            // even winding number -> outside -> sdf is positive
            g.row(v) = pow(-1, w) * g_curr;
        }
        Eigen::MatrixXd V_new = V_curr - h * g;
        Vf_flow.emplace_back(V_new);
        V_curr = V_new;
        iteration++;
    }
    return iteration;
}
