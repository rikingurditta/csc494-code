#ifndef CAGES_BAPN_H
#define CAGES_BAPN_H

#include <Eigen/Dense>

// TODO: use TinyAD for these

// dhat = barrier function safe distance - only need to evaluate distance for pairs that are closer together than dhat
void barrier_aware_projected_newton(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E,
                                    double dhat);

// V, E = input mesh
// dhat = barrier function safe distance
// Chat = constraint set, i.e. list of primitives (edges) that are closer together than dhat
// probably need to use some sort of bounding volume hierarchy for this
void constraint_set(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, double dhat,
                    Eigen::MatrixXi &Chat);


// barrier energy = total energy of system + barrier energy
void barrier_augmented_potential(double dhat, double stiffness);


// V, E = mesh
// x_query = vertex positions for which we want to measure energy
// x = current vertex positions
// v = current velocity
// fe = external forces
// M = mass matrix
// does not take into account friction because we don't need it for nested cages
void total_energy(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E,
                  const Eigen::MatrixXd &x_query,
                  const Eigen::MatrixXd &x, const Eigen::MatrixXd &v, const Eigen::MatrixXd &fe,
                  const Eigen::MatrixXd &M);


// nested cages mesh quality energy
// for now, will just use area energy
// note for implementation - can use green's theorem to integrate area over boundary, i.e. int_S 1 dA = int_∂S x dy
// V, E = mesh
//        (V unnecessary for area energy)
// x_query = current vertex positions
void E_nested_cages(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E,
                    const Eigen::MatrixXd &x_query);

#endif //CAGES_BAPN_H
