#include "flow.h"
#include "grad_sdf.h"
#include "igl/adjacency_list.h"
#include "closest_points.h"

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
    for (int i = 0; i <= num_points; i++) {
        double t = i * 1. / num_points;
        Eigen::RowVector2d g0;
        grad_unsigned_distance(V, E, (1. - t) * p + t * p_e0, g0);
        double l0 = (p - p_e0).norm();
        Eigen::RowVector2d g1;
        grad_unsigned_distance(V, E, (1. - t) * p + t * p_e1, g1);
        double l1 = (p - p_e1).norm();
        g += pow(1. - t, 2) * (g0 / l0 + g1 / l1);
    }
    g.normalize();
}


int query_winding_number(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E,
                         const Eigen::RowVector2d &p) {
    // TODO: not sure if robust implementation
    Eigen::Matrix2d A;
    A.col(0) = Eigen::Vector2d(1, 0);
    int winding = 0;
    for (int e = 0; e < E.rows(); e++) {
        int vc0 = E(e, 0), vc1 = E(e, 1);
        A.col(1) = (V.row(vc1) - V.row(vc0)).transpose();
        Eigen::Vector2d c = p - V.row(vc0);
        if (A.determinant() != 0) {
            Eigen::Vector2d st = A.inverse() * c;
            if (st(0) < 0. and 0. < st(1) and st(1) < 1.)
                winding++;
        }
    }
    return winding;
}

bool query_point_inside(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::Vector2d &p) {
    return query_winding_number(V, E, p) % 2 == 1;
}

bool query_mesh_inside(const Eigen::MatrixXd &Vc, const Eigen::MatrixXi &Ec,
                       const Eigen::MatrixXd &Vf, const Eigen::MatrixXi &Ef) {
    // TODO: this only checks if all vertices are inside, not if entire mesh is inside
    for (int v = 0; v < Vf.rows(); v++) {
        if (!query_point_inside(Vc, Ec, Vf.row(v)))
            return false;
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
    while (!query_mesh_inside(Vc, Ec, Vf, Ef) && iteration < max_num_meshes) {
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
