#ifndef CAGES_CLOSEST_POINTS_H
#define CAGES_CLOSEST_POINTS_H

#include <Eigen/Dense>

void closest_point_on_segment(const Eigen::Vector2d &e0, const Eigen::Vector2d &e1, const Eigen::Vector2d &p,
                              Eigen::Vector2d &closest_point, double &d_sq);

void closest_point_on_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::RowVector2d &p,
                           Eigen::RowVector2d &closest_point);

void query_closest_points(const Eigen::MatrixXd &Vf, const Eigen::MatrixXd &Vc, const Eigen::MatrixXi &Ec,
                          Eigen::MatrixXd &C);

#endif //CAGES_CLOSEST_POINTS_H
