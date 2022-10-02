#ifndef CAGES_FLOW_H
#define CAGES_FLOW_H

#include "Eigen/Dense"

void flow(Eigen::MatrixXd &Vc, Eigen::MatrixXi &Ec,
          Eigen::MatrixXd &Vf, Eigen::MatrixXi &Ef,
          std::vector<Eigen::MatrixXd> &Vf_flow);

#endif //CAGES_FLOW_H
