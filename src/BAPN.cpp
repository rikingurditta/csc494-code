#include <vector>
#include <tuple>
#include "BAPN.h"
#include <TinyAD/Utils/Helpers.hh>
#include <TinyAD/ScalarFunction.hh>

double clamp(double val, double lo, double hi) {
    return std::max(std::min(val, lo), hi);
}

void barrier_aware_projected_newton(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E,
                                    double dhat);

// V0, E0, V1, E1 = current meshes
// dhat = barrier function safe distance
// Chat = constraint set, i.e. list of primitives (edges) that are closer together than dhat
// probably need to use some sort of bounding volume hierarchy for this
std::vector<std::tuple<int, int, int>> constraint_set(const Eigen::MatrixXd &V0, const Eigen::MatrixXi &E0,
                                                      const Eigen::MatrixXd &V1,
                                                      double dhat) {
    // TODO: redo with BVH
    auto out = std::vector<std::tuple<int, int, int>>();
    for (int e = 0; e < E0.rows(); e++) {
        int ve0 = E0(e, 0), ve1 = E0(e, 1);
        Eigen::Vector2d e_dir = V0.row(ve1) - V0.row(ve0);
        Eigen::Matrix2d M;
        M.col(0) = e_dir;
        M.col(1) << e_dir(1), -e_dir(0);
        if (M.determinant() == 0)
            continue;
        for (int v = 0; v < V1.rows(); v++) {
            Eigen::Vector2d v_curr = V1.row(v);
            Eigen::Vector2d st = M.inverse() * v_curr;
            st(0) = clamp(st(0), 0, 1);
            double sqdist = (v_curr - M * st).squaredNorm();
            if (sqdist < dhat * dhat) {
                out.emplace_back(ve0, ve1, v);
            }
        }
    }
    return out;
}

std::tuple<double, Eigen::VectorXd, Eigen::SparseMatrix<double>> nested_cages_energy(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::VectorXd &x_query) {
    auto func = TinyAD::scalar_function<2>(TinyAD::range(x_query.rows() / 2));

    // put together 2d vertex positions from global x vector
    Eigen::VectorXd x = func.x_from_data([&](int v) {
        return Eigen::Vector2d(x_query(v * 2), x_query(v * 2 + 1));
    });
    func.add_elements<2>(TinyAD::range(E.rows()), [&](auto &element) {
        using T = TINYAD_SCALAR_TYPE(element);

        Eigen::Index e = element.handle;
        std::cout << e << "\n";
        std::cout << element.variables(E(e, 0)) << "\n";
        Eigen::Vector2<T> v0 = element.variables(E(e, 0));
        Eigen::Vector2<T> v1 = element.variables(E(e, 1));
        return 0.5 * (v1.y() - v0.y()) * (v0.x() + v1.x());
    });
    return func.eval_with_hessian_proj(x_query);
}

auto barrier_potential(const Eigen::MatrixXd &V0, const Eigen::MatrixXi &E0,
                       const Eigen::MatrixXd &V1, const Eigen::MatrixXi &E1,
                       const Eigen::VectorXd &x_query,
                       double dhat, double stiffness) {
    auto Chat = constraint_set(V0, E0, V1, dhat);
    auto func = TinyAD::scalar_function<2>(TinyAD::range(x_query.rows() / 2));

    func.x_from_data([&](int v) {
        return Eigen::Vector2d(x_query(v * 2), x_query(v * 2 + 1));
    });

    func.add_elements<3>(Chat, [&](auto &element) {
        using T = TINYAD_SCALAR_TYPE(element);
        auto [i_e0, i_e1, i_v] = element.handle;
        Eigen::Vector2<T> v = element.variables(i_v);
        Eigen::Vector2<T> e_dir = element.variables(i_e1) - element.variables(i_e0);
        Eigen::Matrix2<T> M;
        M.col(0) = e_dir;
        M.col(1) << e_dir(1), -e_dir(0);
        if (M.determinant() == 0)
            return (T) 0;
        Eigen::Vector2<T> st = M.inverse() * v;
        st(0) = clamp(st(0), (T) 0, (T) 1);
        T dist = (v - M * st).norm();
        return -(dist - dhat) * log(dist / dhat);
    });

    return func.eval_with_hessian_proj(x_query);
}

auto total_energy(const Eigen::MatrixXd &Vc, const Eigen::MatrixXi &Ec,
                  const Eigen::MatrixXd &Vf, const Eigen::MatrixXi &Ef,
                  const Eigen::VectorXd &x, const Eigen::VectorXd &x_query,
                  const Eigen::VectorXd &v, const Eigen::VectorXd &fe,
                  const Eigen::MatrixXd &M, double dt) {
    double dhat = 0.1;
    double lambda = 0.1;
    VectorXAD xx = x_query;
//    Eigen::VectorXd xhat = x + dt * v + dt * dt * M.inverse() * fe;
//    ADouble E = (xx - xhat).transpose() * M * (xx - xhat);
    auto EgH_cages = nested_cages_energy(Vc, Ec, x_query);
    auto EgH_barrier_f = barrier_potential(Vc, Ec, Vf, Ef, x_query, dhat, 0.1);
    auto EgH_barrier_c = barrier_potential(Vf, Ef, Vc, Ef, x_query, dhat, 0.1);

    VectorXAD xx_f = xx.segment(Vc.rows() * 2, Vf.rows() * 2);
    ADouble E_disp = (xx_f - x.segment(Vc.rows() * 2, Vf.rows() * 2)).squaredNorm() / 2. * lambda;

    double energy = std::get<0>(EgH_cages)
                    + std::get<0>(EgH_barrier_f) + std::get<0>(EgH_barrier_c)
                    + E_disp.val;
    Eigen::VectorXd gradient = std::get<1>(EgH_cages)
                               + std::get<1>(EgH_barrier_f) + std::get<1>(EgH_barrier_c)
                               + E_disp.grad;
    Eigen::SparseMatrix<double> hessian_proj = std::get<2>(EgH_cages)
                                               + std::get<2>(EgH_barrier_f) + std::get<2>(EgH_barrier_c)
                                               // don't need to proj E_disp.Hess since its hessian is already SPD
                                               + E_disp.Hess;
    return std::tuple<double, Eigen::VectorXd, Eigen::SparseMatrix<double>>(energy, gradient, hessian_proj);
}
