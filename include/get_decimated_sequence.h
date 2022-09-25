#ifndef CAGES_GET_DECIMATED_SEQUENCE_H
#define CAGES_GET_DECIMATED_SEQUENCE_H

#include "Eigen/Dense"
#include "igl/decimate.h"

bool get_decimated_sequence(Eigen::MatrixXd &V, Eigen::MatrixXi &F, int num_meshes, int coarsest_faces,
                            std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> &out);

#endif //CAGES_GET_DECIMATED_SEQUENCE_H
