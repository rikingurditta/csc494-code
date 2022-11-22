#include <vector>
#include <tuple>
#include "BAPN.h"
#include <TinyAD/Utils/Helpers.hh>
#include <TinyAD/ScalarFunction.hh>
#include "helpers.h"


void barrier_aware_projected_newton(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E,
                                    double dhat);

// V0, E0, V1, E1 = current meshes
// dhat = barrier function safe distance
// Chat = constraint set, i.e. list of primitives (edges) that are closer together than dhat
// probably need to use some sort of bounding volume hierarchy for this
std::vector<std::tuple<int, int, int>> constraint_set(const Eigen::MatrixXi &E, const Eigen::VectorXd &x_edges,
                                                      const Eigen::VectorXd &x_vertices, double dhat) {
    // TODO: redo with BVH
    auto out = std::vector<std::tuple<int, int, int>>();
    for (int e = 0; e < E.rows(); e++) {
        int e0 = E(e, 0), e1 = E(e, 1);
        Eigen::Vector2d ve0 = x_edges.segment<2>(e0 * 2);
        Eigen::Vector2d ve1 = x_edges.segment<2>(e1 * 2);
        Eigen::Vector2d e_dir = ve1 - ve0;
        Eigen::Matrix2d M;
        M.col(0) = e_dir;
        M.col(1) << e_dir(1), -e_dir(0);
        if (M.determinant() == 0)
            continue;
        for (int v = 0; v < x_vertices.rows() / 2; v++) {
            Eigen::Vector2d v_curr = x_vertices.segment<2>(v * 2);
            Eigen::Vector2d st = M.inverse() * (v_curr - ve0);
            double s = std::clamp(st(0), 0., 1.);
            double sqdist = (v_curr - (ve0 + s * e_dir)).squaredNorm();
            if (sqdist < dhat * dhat) {
                out.emplace_back(e0, e1, v);
            }
        }
    }
    return out;
}

std::tuple<double, Eigen::VectorXd, Eigen::SparseMatrix<double>>
nested_cages_energy(const Eigen::MatrixXi &E, const Eigen::VectorXd &x_query, int index_offset) {
    auto func = TinyAD::scalar_function<2>(TinyAD::range(x_query.rows() / 2 - index_offset));

    // put together 2d vertex positions from global x vector
    Eigen::VectorXd x = func.x_from_data([&](int v) {
        return Eigen::Vector2d(x_query((v + index_offset) * 2), x_query((v + index_offset) * 2 + 1));
    });
    func.add_elements<2>(TinyAD::range(E.rows()), [&](auto &element) {
        using T = TINYAD_SCALAR_TYPE(element);

        Eigen::Index e = element.handle;
        Eigen::Vector2<T> v0 = element.variables(E(e, 0));
        Eigen::Vector2<T> v1 = element.variables(E(e, 1));
        return -0.5 * (v1.y() - v0.y()) * (v0.x() + v1.x());
    });
    return func.eval_with_hessian_proj(x_query);
}

// assumes x = stack(flatten(V0), flatten(V1))
// typically, x = stack(flatten(Vc), flatten(Vf)), so if want to use Vc for vertices, set use_mesh0_edges = false
std::tuple<double, Eigen::VectorXd, Eigen::SparseMatrix<double>>
barrier_potential(const Eigen::MatrixXi &E, const Eigen::VectorXd &x_edges,
                  const Eigen::VectorXd &x_vertices,
                  double dhat, double stiffness) {
    std::vector<std::tuple<int, int, int>> Chat;
    Chat = constraint_set(E, x_edges, x_vertices, dhat);

    auto func = TinyAD::scalar_function<2>(TinyAD::range(x_vertices.rows() / 2));

    func.x_from_data([&](int v) {
        return Eigen::Vector2d(x_vertices(v * 2), x_vertices(v * 2 + 1));
    });

    func.add_elements<1>(Chat, [&](auto &element) {
        using T = TINYAD_SCALAR_TYPE(element);
        auto [i_e0, i_e1, i_v] = element.handle;
        Eigen::Vector2<T> v = element.variables(i_v);
        Eigen::Vector2d ve0 = x_edges.segment<2>(i_e0 * 2);
        Eigen::Vector2d ve1 = x_edges.segment<2>(i_e1 * 2);
        Eigen::Vector2d e_dir = ve1 - ve0;
        Eigen::Matrix2d M;
        M.col(0) = e_dir;
        M.col(1) << e_dir(1), -e_dir(0);
        if (M.determinant() == 0)
            return (T) 0;
        Eigen::Vector2<T> st = M.inverse() * (v - ve0);
        T s = std::clamp(st(0), (T) 0., (T) 1.);
        T dist = (v - (ve0 + s * e_dir)).norm();
        if (abs(dist) < 1e-15)
            return (T) 0;
        else
            return -(dist - dhat) * log(dist / dhat);
    });

    return func.eval_with_hessian_proj(x_vertices);
}

// x = stack(flatten(Vc), flatten(Vf))
std::tuple<double, Eigen::VectorXd, Eigen::SparseMatrix<double>>
total_energy(const Eigen::MatrixXd &Vc, const Eigen::MatrixXi &Ec, const Eigen::MatrixXd &Vf,
             const Eigen::MatrixXd &Vf_next, const Eigen::MatrixXi &Ef, const Eigen::VectorXd &x) {
    double dhat = 0.1;
    double lambda = 10;
    double k = 0.1;

    Eigen::VectorXd xc = x.segment(0, Vc.rows() * 2);
    Eigen::VectorXd xf = x.segment(Vc.rows() * 2, Vf.rows() * 2);
    Eigen::VectorXd xf_next = flatten(Vf_next);
    Eigen::Index xc_start = 0;
    Eigen::Index xf_start = xc.size();

    auto [E_cages, g_cages, H_cages]
            = nested_cages_energy(Ec, x, 0);
    // compare vertices of F with edges of C
//    auto [E_barrier_F_vs_C, gf_barrier_F_vs_C, Hf_barrier_F_vs_C]
//            = barrier_potential(Ec, xc, xf, dhat, k);
//    Eigen::VectorXd g_barrier_F_vs_C = vector_segment_assemble(x.size(), xf_start, gf_barrier_F_vs_C);
//    Eigen::SparseMatrix<double> H_barrier_F_vs_C =
//            sparse_block_assemble(x.size(), x.size(), xf_start, xf_start, Hf_barrier_F_vs_C);
    // compare vertices of C with edges of F
//    auto [E_barrier_C_vs_F, gc_barrier_C_vs_F, Hc_barrier_C_vs_F]
//            = barrier_potential(Ef, xf, xc, dhat, k);
//    Eigen::VectorXd g_barrier_C_vs_F = vector_segment_assemble(x.size(), xc_start, gc_barrier_C_vs_F);
//    Eigen::SparseMatrix<double> H_barrier_C_vs_F =
//            sparse_block_assemble(x.size(), x.size(), xc_start, xc_start, Hc_barrier_C_vs_F);
    // compare vertices of C with edges of C
    auto [E_barrier_C_vs_C, gc_barrier_C_vs_C, Hc_barrier_C_vs_C]
            = barrier_potential(Ec, xc, xc, dhat, k);
    Eigen::VectorXd g_barrier_C_vs_C = vector_segment_assemble(x.size(), xc_start, gc_barrier_C_vs_C);
    Eigen::SparseMatrix<double> H_barrier_C_vs_C =
            sparse_block_assemble(x.size(), x.size(), xc_start, xc_start, Hc_barrier_C_vs_C);

    // calculating displacement energy and its derivatives manually because TinyAD is annoying when it comes to dynamic
    // size matrices
    double E_disp = (xf - xf_next).squaredNorm() / 2. * lambda;
    Eigen::VectorXd g_disp_f = (xf - xf_next) * lambda;
    Eigen::VectorXd g_disp = Eigen::VectorXd::Zero(x.rows());
    g_disp.segment(Vc.rows() * 2, Vf.rows() * 2) = g_disp_f;
    // don't need to proj E_disp.Hess since its hessian is already SPD
    Eigen::MatrixXd H_disp_f = Eigen::MatrixXd::Identity(xf.rows(), xf.rows()) * lambda;
    Eigen::MatrixXd H_disp = Eigen::MatrixXd::Zero(x.rows(), x.rows());
    H_disp.block(Vc.rows() * 2, Vc.rows() * 2, Vf.rows() * 2, Vf.rows() * 2) = H_disp_f;

//    double energy = E_cages + E_barrier_F_vs_C + E_barrier_C_vs_F + E_barrier_C_vs_C + E_disp;
//    Eigen::VectorXd gradient = g_cages + g_barrier_F_vs_C + g_barrier_C_vs_F + g_barrier_C_vs_C + g_disp;
//    Eigen::SparseMatrix<double> hessian_proj =
//            H_cages + H_barrier_F_vs_C + H_barrier_C_vs_F + H_barrier_C_vs_C + H_disp;

//    return {energy, gradient, hessian_proj};
//    return {E_cages, g_cages, H_cages};
    return {E_cages + E_barrier_C_vs_C, g_cages + g_barrier_C_vs_C, H_cages + H_barrier_C_vs_C};
}
