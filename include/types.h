#ifndef CAGES_EIGENTYPES_H
#define CAGES_EIGENTYPES_H

#include <Scalar.hh>

using ADouble = TinyAD::Double<1>;
using MatrixXAD = Eigen::Matrix<ADouble, Eigen::Dynamic, Eigen::Dynamic>;

#endif //CAGES_EIGENTYPES_H
