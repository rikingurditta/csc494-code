#ifndef CAGES_BAPN_H
#define CAGES_BAPN_H

#include <Eigen/Dense>

// TODO: use TinyAD for these

// V0, E0, V1, E1 = current meshes
// dhat = barrier function safe distance - only need to evaluate distance for pairs that are closer together than dhat
void barrier_aware_projected_newton(const Eigen::MatrixXd &V0, const Eigen::MatrixXi &E0,
                                    const Eigen::MatrixXd &V1, const Eigen::MatrixXi &E1,
                                    double dhat);

// V0, E0, V1, E1 = current meshes
// dhat = barrier function safe distance
// Chat = constraint set, i.e. list of primitives (edges) that are closer together than dhat
// probably need to use some sort of bounding volume hierarchy for this
void constraint_set(const Eigen::MatrixXd &V0, const Eigen::MatrixXi &E0,
                    const Eigen::MatrixXd &V1, const Eigen::MatrixXi &E1,
                    double dhat,
                    Eigen::MatrixXi &Chat);


// V, E = current meshe
// x = current mesh vertex positions
// x_query = vertex positions for which the energy is being calculated
// v = current mesh velocity
// fe = external forces
// M = mass matrix
// does not take into account friction because we don't need it for nested cages
// returns nested cages energy + barrier energy
double total_energy(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E,
                    const Eigen::VectorXd &x, const Eigen::VectorXd &x_query,
                    const Eigen::VectorXd &v, const Eigen::VectorXd &fe,
                    const Eigen::MatrixXd &M,
                    double dt);


// barrier energy = total energy of system + barrier energy
// V0, E0, V1, E1 = current meshes
double barrier_potential(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E,
                         const Eigen::VectorXd &x_query,
                         double dhat, double stiffness);


// nested cages mesh quality energy
// for now, will just use area energy
// note for implementation - can use green's theorem to integrate area over boundary, i.e. int_S 1 dA = int_∂S x dy
// V, E = current mesh
// x_query = vertex positions for which the energy is being calculated
// returns energy meant to be minimized for nested cages (e.g. volume energy in a typical case)
double nested_cages_energy(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::VectorXd &x_query);

#endif //CAGES_BAPN_H
