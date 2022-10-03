#include <utility>
#include <iostream>
#include "flow.h"
#include "grad_sdf.h"

int query_winding_number(Eigen::MatrixXd &V, Eigen::MatrixXi &E, const Eigen::RowVector2d &p)
{
    // TODO: not sure if robust implementation
    Eigen::Matrix2d A;
    A.col(0) = Eigen::Vector2d(1, 0);
    int winding = 0;
    for (int e = 0; e < E.rows(); e++)
    {
        int vc0 = E(e, 0), vc1 = E(e, 1);
        A.col(1) = (V.row(vc1) - V.row(vc0)).transpose();
        Eigen::Vector2d c = p - V.row(vc0);
        if (A.determinant() != 0)
        {
            Eigen::Vector2d st = A.inverse() * c;
            if (st(0) < 0. and 0. < st(1) and st(1) < 1.)
                winding++;
        }
    }
    return winding;
}

bool query_point_inside(Eigen::MatrixXd &V, Eigen::MatrixXi &E, const Eigen::Vector2d &p)
{
    return query_winding_number(V, E, p) % 2 == 1;
}

bool query_mesh_inside(Eigen::MatrixXd &Vc, Eigen::MatrixXi &Ec,
                       Eigen::MatrixXd &Vf, Eigen::MatrixXi &Ef)
{
    // TODO: this only checks if all vertices are inside, not if entire mesh is inside
    for (int v = 0; v < Vf.rows(); v++)
    {
        if (!query_point_inside(Vc, Ec, Vf.row(v)))
            return false;
    }
    return true;
}

void flow(Eigen::MatrixXd &Vc, Eigen::MatrixXi &Ec,
          Eigen::MatrixXd &Vf, Eigen::MatrixXi &Ef,
          std::vector<Eigen::MatrixXd> &Vf_flow)
{
    int c = 0;
    double h = 0.1;
    Vf_flow.emplace_back(Vf);
    Eigen::MatrixXd V_curr = Vf;
    while (!query_mesh_inside(Vc, Ec, Vf, Ef) && c < 25) {
        // g = ∇ sdf
        // TODO: calculate SDF better, use numerical integration along edges to smooth it
        Eigen::MatrixXd g = Eigen::MatrixXd::Zero(Vf.rows(), Vf.cols());
        for (int v = 0; v < V_curr.rows(); v++) {
            int min_index = 0;
            double min_sqdistance = (V_curr.row(v) - Vc.row(0)).squaredNorm();
            for (int vc = 1; vc < Vc.rows(); vc++) {
                double sqdistance = (V_curr.row(v) - Vc.row(vc)).squaredNorm();
                if (sqdistance < min_sqdistance) {
                    min_index = vc;
                    min_sqdistance = sqdistance;
                }
            }
            // calculate winding number
            int w = query_winding_number(Vc, Ec, V_curr.row(v));
            // even winding number -> outside -> sdf is positive
            g.row(v) = pow(-1, w)
                       * (V_curr.row(v) - Vc.row(min_index)).normalized();
        }
        Eigen::MatrixXd V_new = V_curr - h * g;
        Vf_flow.emplace_back(V_new);
        V_curr = V_new;
        std::cout << c << "\n";
        c++;
    }
}
