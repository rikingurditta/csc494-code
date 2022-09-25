#ifndef CAGES_FLOW_H
#define CAGES_FLOW_H

#include "Eigen/Dense"

void flow(Eigen::MatrixXd &Vc, Eigen::MatrixXi &Fc,
          Eigen::MatrixXd &Vf, Eigen::MatrixXd &Ff,
          std::vector<Eigen::MatrixXd> &Vf_flow);

#endif //CAGES_FLOW_H
