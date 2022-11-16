#ifndef CAGES_EIGENTYPES_H
#define CAGES_EIGENTYPES_H

#include <Eigen/Dense>
#include <TinyAD/Scalar.hh>

using ADouble = TinyAD::Double<1>;
using VectorXAD = Eigen::Vector<ADouble, Eigen::Dynamic>;
using MatrixXAD = Eigen::Matrix<ADouble, Eigen::Dynamic, Eigen::Dynamic>;

#endif //CAGES_EIGENTYPES_H
