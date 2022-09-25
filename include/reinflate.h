#ifndef CAGES_REINFLATE_H
#define CAGES_REINFLATE_H

#include "Eigen/Dense"

void reinflate(Eigen::MatrixXd Vc, Eigen::MatrixXi Fc,
               Eigen::MatrixXd Vf, Eigen::MatrixXi Ff,
               std::vector<Eigen::MatrixXd> &Vf_flow,
               Eigen::MatrixXd Vc_inflated);

void reinflate_simulate_timestep(Eigen::MatrixXd Vc, Eigen::MatrixXi Fc,
                                 Eigen::MatrixXd Vf, Eigen::MatrixXi Ff,
                                 double t,
                                 Eigen::MatrixXd Vc_inflated);

#endif //CAGES_REINFLATE_H