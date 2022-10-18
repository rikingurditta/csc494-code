#ifndef CAGES_REINFLATE_H
#define CAGES_REINFLATE_H

#include "Eigen/Dense"
#include <vector>

void reinflate(const Eigen::MatrixXd &Vc, const Eigen::MatrixXi &Ec,
               const Eigen::MatrixXd &Vf, const Eigen::MatrixXi &Ef,
               const std::vector<Eigen::MatrixXd> &Vf_flow,
               Eigen::MatrixXd &Vc_inflated);

void reinflate_simulate_timestep(const Eigen::MatrixXd &Vc, const Eigen::MatrixXi &Ec,
                                 const Eigen::MatrixXd &Vf_curr, const Eigen::MatrixXd &Vf_next, const Eigen::MatrixXi &Ef,
                                 Eigen::MatrixXd &Vc_inflated);

#endif //CAGES_REINFLATE_H
