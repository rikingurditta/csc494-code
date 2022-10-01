#ifndef CAGES_SDF_FIELD_H
#define CAGES_SDF_FIELD_H

#include "Eigen/Dense"

void grad_sdf(Eigen::MatrixXd &Vc, Eigen::MatrixXi &Ec,
              Eigen::MatrixXd &Vf, Eigen::MatrixXd &Ef,
              Eigen::MatrixXd &g);

#endif //CAGES_SDF_FIELD_H
