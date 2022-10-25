#include <vector>
#include "BAPN.h"

double clamp(double val, double lo, double hi)
{
    return std::max(std::min(val, lo), hi);
}

void barrier_aware_projected_newton(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E,
                                    double dhat);

// V0, E0, V1, E1 = current meshes
// dhat = barrier function safe distance
// Chat = constraint set, i.e. list of primitives (edges) that are closer together than dhat
// probably need to use some sort of bounding volume hierarchy for this
void constraint_set(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E,
                    double dhat,
                    Eigen::MatrixXi &Chat) {
    Chat = Eigen::MatrixXi::Zero(E.rows() * (E.rows() + 1) / 2, 2);
    for (int e0 = 0; e0 < E.rows(); e0++) {
        for (int e1 = e0 + 1; e1 < E.rows(); e1++) {
            Chat(e0 * E.rows() + e1, 0) = e0;
            Chat(e0 * E.rows() + e1, 1) = e1;
        }
    }
}


double total_energy(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E,
                    const Eigen::VectorXd &x, const Eigen::VectorXd &x_query,
                    const Eigen::VectorXd &v, const Eigen::VectorXd &fe,
                    const Eigen::MatrixXd &M,
                    double dt) {
    Eigen::VectorXd xhat = x + dt * v + dt * dt * M.inverse() * fe;
    return 0.5 * (x_query - xhat).transpose() * M * (x_query - xhat)
           + dt * dt * nested_cages_energy(V, E, x_query)
           + barrier_potential(V, E, x_query, 0.1, 0.1);
}


double barrier_potential(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E,
                         const Eigen::VectorXd &x_query,
                         double dhat, double stiffness) {
    double energy = 0;
    Eigen::MatrixXi Chat;
    constraint_set(V, E, dhat, Chat);
    for (int pair = 0; pair < Chat.rows(); pair++) {
        int e0 = Chat(pair, 0);
        int e1 = Chat(pair, 1);
        Eigen::Vector2d a0 = x_query.segment(E(e0, 0) * 2, 2);
        Eigen::Vector2d b0 = x_query.segment(E(e0, 1) * 2, 2);
        Eigen::Vector2d a1 = x_query.segment(E(e1, 0) * 2, 2);
        Eigen::Vector2d b1 = x_query.segment(E(e1, 1) * 2, 2);
        Eigen::Vector2d d0 = b0 - a0;
        Eigen::Vector2d d1 = b1 - a1;
        double a = d0.x(), b = d0.y(), c = d1.x(), d = d1.y();
        Eigen::Matrix2d A;
        A << a, c, b, d;
        Eigen::Vector2d B;
        B << (a1 - a0).dot(d0), (a1 - a0).dot(d1);
        if (A.determinant() != 0)
        {
            Eigen::Vector2d st = A.inverse() * B;
            double s = clamp(st(0), 0, 1);
            double t = clamp(st(1), 0, 1);
            double dist = (a0 + s * d0 - a1 - t * d1).norm();
            energy += -(dist - dhat) * log(dist / dhat);
        }
    }
}


double nested_cages_energy(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::VectorXd &x_query) {
    // assuming shape is a polygon with no holes (boundary is connected) and all edges are oriented correctly
    // TODO: implement weird topology e.g. holes using winding number

    // using green's theorem for the special case of piecewise linear boundary,
    // area is sum over edges of 0.5 * (y1 - y0) * (x0 + x1)
    double area = 0;
    for (int e = 0; e < E.rows(); e++) {
        int v0 = E(e, 0), v1 = E(e, 1);
        double x0 = x_query(v0 * 2), y0 = x_query(v0 * 2 + 1);
        double x1 = x_query(v1 * 2), y1 = x_query(v1 * 2 + 1);
        area += 0.5 * (y1 - y0) * (x0 + x1);
    }
}
