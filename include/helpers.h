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

#endif //CAGES_HELPERS_H
