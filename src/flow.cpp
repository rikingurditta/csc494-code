#include "flow.h"

void flow(Eigen::MatrixXd &Vc, Eigen::MatrixXi &Fc,
          Eigen::MatrixXd &Vf, Eigen::MatrixXd &Ff,
          std::vector<Eigen::MatrixXd> &Vf_flow)
{
    Eigen::MatrixXd g = Eigen::MatrixXd::Zero(Vf.rows(), Vf.cols());
}
