#include <iostream>
#include "flow.h"
#include "grad_sdf.h"

int winding(Eigen::MatrixXd V, Eigen::MatrixXi E, Eigen::RowVector2d p)
{
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

void flow(Eigen::MatrixXd &Vc, Eigen::MatrixXi &Ec,
          Eigen::MatrixXd &Vf, Eigen::MatrixXi &Ef,
          std::vector<Eigen::MatrixXd> &Vf_flow)
{
    Eigen::MatrixXd g = Eigen::MatrixXd::Zero(Vf.rows(), Vf.cols());
    for (int vf = 0; vf < Vf.rows(); vf++)
    {
        int min_index = 0;
        double min_sqdistance = (Vf.row(vf) - Vc.row(0)).squaredNorm();
        for (int vc = 1; vc < Vc.rows(); vc++)
        {
            double sqdistance = (Vf.row(vf) - Vc.row(vc)).squaredNorm();
            if (sqdistance < min_sqdistance)
            {
                min_index = vc;
                min_sqdistance = sqdistance;
            }
        }
        // calculate winding number
        int w = winding(Vc, Ec, Vf.row(vf));
        bool inside = w % 2 == 1;
        std::cout << "vf: " << Vf.row(vf) << " winding: " << w << "\n";
    }
}
