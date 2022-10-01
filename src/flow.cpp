#include "flow.h"

void flow(Eigen::MatrixXd &Vc, Eigen::MatrixXi &Ec,
          Eigen::MatrixXd &Vf, Eigen::MatrixXd &Ef,
          std::vector<Eigen::MatrixXd> &Vf_flow)
{
    Eigen::MatrixXd g = Eigen::MatrixXd::Zero(Vf.rows(), Vf.cols());
}
