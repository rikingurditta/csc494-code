#ifndef CAGES_HELPERS_H
#define CAGES_HELPERS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

Eigen::VectorXd flatten(const Eigen::MatrixXd &V);

Eigen::MatrixXd unflatten(const Eigen::VectorXd &v, int dim);

Eigen::VectorXd vector_segment_assemble(Eigen::Index size, Eigen::Index start, const Eigen::VectorXd &segment);

Eigen::SparseMatrix<double> sparse_block_assemble(Eigen::Index rows, Eigen::Index cols,
                                                  Eigen::Index startRow, Eigen::Index startCol,
                                                  const Eigen::SparseMatrix<double> &block);

int query_winding_number(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E,
                         const Eigen::RowVector2d &p);

bool query_point_inside(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::Vector2d &p);

bool query_mesh_inside(const Eigen::MatrixXd &V_outside, const Eigen::MatrixXi &E_outside,
                       const Eigen::MatrixXd &V_inside, const Eigen::MatrixXi &E_inside);

bool query_meshes_intersect(const Eigen::VectorXd &x0, const Eigen::MatrixXi &E0,
                            const Eigen::VectorXd &x1, const Eigen::MatrixXi &E1);

#endif //CAGES_HELPERS_H
