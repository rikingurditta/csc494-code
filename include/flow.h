#ifndef CAGES_FLOW_H
#define CAGES_FLOW_H

#include <vector>
#include "Eigen/Dense"

void closest_point_on_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::RowVector2d &p,
                           Eigen::RowVector2d &closest_point);

int flow(const Eigen::MatrixXd &Vc, const Eigen::MatrixXi &Ec, const Eigen::MatrixXd &Vf, const Eigen::MatrixXi &Ef,
         int max_flow_steps, std::vector<Eigen::MatrixXd> &Vf_flow, double h);

#endif //CAGES_FLOW_H
