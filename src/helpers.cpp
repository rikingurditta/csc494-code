#include "helpers.h"

Eigen::VectorXd flatten(Eigen::MatrixXd &V)
{
    Eigen::VectorXd out = Eigen::VectorXd::Zero(V.rows() * V.cols());
    for (int r = 0; r < V.rows(); r++)
    {
        out.segment(r * V.cols(), V.cols()) = V.row(r);
    }
    return out;
}
