#ifndef CAGES_BAPN_H
#define CAGES_BAPN_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

// TODO: use TinyAD for these

// V0, E0, V1, E1 = current meshes
// dhat = barrier function safe distance - only need to evaluate distance for pairs that are closer together than dhat
void barrier_aware_projected_newton(const Eigen::MatrixXd &V0, const Eigen::MatrixXi &E0,
                                    const Eigen::MatrixXd &V1, const Eigen::MatrixXi &E1,
                                    double dhat);

// V0, E0, V1, E1 = current meshes
// dhat = barrier function safe distance
// return constraint set, i.e. list of primitives (edges) that are closer together than dhat
std::vector<std::tuple<int, int, int>> constraint_set(const Eigen::MatrixXi &E, const Eigen::VectorXd &x_edges,
                                                      const Eigen::VectorXd &x_vertices, double dhat);


// V, E = current meshes
// x = queried mesh vertex positions
// does not take into account friction because we don't need it for nested cages
// returns nested cages energy + barrier nested_cages_energy
std::tuple<double, Eigen::VectorXd, Eigen::SparseMatrix<double>>
total_energy(const Eigen::MatrixXd &Vc, const Eigen::MatrixXi &Ec, const Eigen::MatrixXd &Vf,
             const Eigen::MatrixXd &Vf_next, const Eigen::MatrixXi &Ef, const Eigen::VectorXd &x, double dhat);


std::tuple<double, Eigen::VectorXd, Eigen::SparseMatrix<double>>
barrier_potential(const Eigen::MatrixXi &E, const Eigen::VectorXd &x,
                  Eigen::Index e_mesh_offset, Eigen::Index e_mesh_n, Eigen::Index v_mesh_offset, Eigen::Index v_mesh_n,
                  double dhat);


// nested cages mesh quality nested_cages_energy
// for now, will just use area nested_cages_energy
// note for implementation - can use green's theorem to integrate area over boundary, i.e. int_S 1 dA = int_∂S x dy
// V, E = current mesh
// x_query = vertex positions for which the nested_cages_energy is being calculated
// returns energy meant to be minimized for nested cages (e.g. volume nested_cages_energy in a typical case)
std::tuple<double, Eigen::VectorXd, Eigen::SparseMatrix<double>>
nested_cages_energy(const Eigen::MatrixXi &E, const Eigen::VectorXd &x_query, int index_offset);

#endif //CAGES_BAPN_H
