#include "closest_points.h"


void closest_point_on_segment(const Eigen::Vector2d &e0, const Eigen::Vector2d &e1, const Eigen::Vector2d &p,
                              Eigen::Vector2d &closest_point, double &d_sq) {
    Eigen::Vector2d edge_dir = e1 - e0;
    Eigen::Vector2d edge_normal(-edge_dir(1), edge_dir(0));
    Eigen::Matrix2d A;
    A.col(0) = edge_dir;
    A.col(1) = edge_normal;
    // det is only 0 when e0 = e1
    if (A.determinant() == 0) {
        closest_point = e0;
    } else {
        Eigen::Vector2d coeffs = A.inverse() * (p - e0);
        double t = coeffs(0);
        if (t < 0) {
            closest_point = e0;
        } else if (t > 1) {
            closest_point = e1;
        } else {
            closest_point = e0 + t * edge_dir;
        }
    }
    d_sq = (p - closest_point).squaredNorm();
}

void closest_point_on_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::RowVector2d &p,
                           Eigen::RowVector2d &closest_point) {
    closest_point = V.row(0);
    double min_sqdist = (p - closest_point).squaredNorm();
    for (int e = 0; e < E.rows(); e++) {
        Eigen::Vector2d c;
        double d_sq;
        closest_point_on_segment(V.row(E(e, 0)), V.row(E(e, 1)), p,
                                 c, d_sq);
        if (d_sq < min_sqdist) {
            closest_point = c;
            min_sqdist = d_sq;
        }
    }
}

void query_closest_points(const Eigen::MatrixXd &Vf, const Eigen::MatrixXd &Vc, const Eigen::MatrixXi &Ec,
                          Eigen::MatrixXd &C) {
    for (int v = 0; v < Vf.rows(); v++) {
        Eigen::RowVector2d c;
        closest_point_on_mesh(Vc, Ec, Vf.row(v), c);
        C.row(v) = c;
    }
}
