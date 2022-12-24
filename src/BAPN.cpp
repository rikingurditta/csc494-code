#include <vector>
#include <tuple>
#include "BAPN.h"
#include <TinyAD/Utils/Helpers.hh>
#include <TinyAD/ScalarFunction.hh>
#include "helpers.h"
#include <limits>


void barrier_aware_projected_newton(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E,
                                    double dhat);

// V0, E0, V1, E1 = current meshes
// dhat = barrier function safe distance
// Chat = constraint set, i.e. list of primitives (edges) that are closer together than dhat
// probably need to use some sort of bounding volume hierarchy for this
std::vector<std::tuple<int, int, int>> constraint_set(const Eigen::MatrixXi &E, const Eigen::VectorXd &x_edges,
                                                      const Eigen::VectorXd &x_vertices, double dhat) {
    // TODO: redo with BVH
    std::vector<std::tuple<int, int, int>> out;
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
    auto func = TinyAD::scalar_function<2>(TinyAD::range(x_query.rows() / 2));

    // put together 2d vertex positions from global x vector
    Eigen::VectorXd x = func.x_from_data([&](int v) {
        return Eigen::Vector2d(x_query((v + index_offset) * 2), x_query((v + index_offset) * 2 + 1));
    });
    func.add_elements<2>(TinyAD::range(E.rows()), [&](auto &element) {
        using T = TINYAD_SCALAR_TYPE(element);

        Eigen::Index e = element.handle;
        Eigen::Vector2<T> v0 = element.variables(E(e, 0));
        Eigen::Vector2<T> v1 = element.variables(E(e, 1));

        // negating because test meshes are cw instead of ccw
        return -0.5 * (v1.y() - v0.y()) * (v0.x() + v1.x());
    });
    auto [Energy, g, H] = func.eval_with_hessian_proj(x_query);
    if (Energy < 0)
        Energy = std::numeric_limits<double>::infinity();
    return {Energy, g, H};
}

std::tuple<double, Eigen::VectorXd, Eigen::SparseMatrix<double>>
barrier_potential(const Eigen::MatrixXi &E, const Eigen::VectorXd &x,
                  Eigen::Index e_mesh_offset, Eigen::Index e_mesh_n, Eigen::Index v_mesh_offset, Eigen::Index v_mesh_n,
                  double dhat) {
    std::vector<std::tuple<int, int, int>> Chat;
    Chat = constraint_set(E,
                          x.segment(e_mesh_offset * 2, e_mesh_n * 2),
                          x.segment(v_mesh_offset * 2, v_mesh_n * 2),
                          dhat);

    auto func = TinyAD::scalar_function<2>(TinyAD::range(x.rows() / 2));

    func.x_from_data([&](int v) {
        return Eigen::Vector2d(x(v * 2), x(v * 2 + 1));
    });

    func.add_elements<3>(Chat, [&](auto &element) {
        using T = TINYAD_SCALAR_TYPE(element);
        auto [i_e0, i_e1, i_v] = element.handle;
        Eigen::Vector2<T> ve0 = element.variables(i_e0 + e_mesh_offset);
        Eigen::Vector2<T> ve1 = element.variables(i_e1 + e_mesh_offset);
        Eigen::Vector2<T> v = element.variables(i_v + v_mesh_offset);
        Eigen::Vector2<T> e_dir = ve1 - ve0;
        Eigen::Matrix2<T> M;
        M.col(0) = e_dir;
        M.col(1) << e_dir(1), -e_dir(0);
        if (M.determinant() == (T) 0)
            return (T) 0;
        Eigen::Vector2<T> st = M.inverse() * (v - ve0);
        T s = std::clamp(st(0), (T) 0., (T) 1.);
        T dist = (v - (ve0 + s * e_dir)).norm();
        if (abs(dist) < 1e-15)
            return (T) 0;
        else
            return -(dist - dhat) * (dist - dhat) * log(dist / dhat);
    });

    return func.eval_with_hessian_proj(x);
}

// x = stack(flatten(Vc), flatten(Vf))
std::tuple<double, Eigen::VectorXd, Eigen::SparseMatrix<double>>
total_energy(const Eigen::MatrixXd &Vc, const Eigen::MatrixXi &Ec, const Eigen::MatrixXd &Vf,
             const Eigen::MatrixXd &Vf_next, const Eigen::MatrixXi &Ef, const Eigen::VectorXd &x, double dhat) {

    auto [E_cages, g_cages, H_cages]
            = nested_cages_energy(Ec, x, 0);

    Eigen::VectorXd xc = x.segment(0, Vc.rows() * 2);
    Eigen::VectorXd xf = x.segment(Vc.rows() * 2, Vf.rows() * 2);
    Eigen::VectorXd xf_next = flatten(Vf_next);
    Eigen::Index xc_start = 0;  // start index for
    Eigen::Index xf_start = xc.size();

    // compare vertices of F with edges of C
    auto [E_barrier_F_vs_C, g_barrier_F_vs_C, H_barrier_F_vs_C] = barrier_potential(Ec, x,
                                                                                    0, Vc.rows(), Vc.rows(), Vf.rows(),
                                                                                    dhat);
    // compare vertices of C with edges of F
    auto [E_barrier_C_vs_F, g_barrier_C_vs_F, H_barrier_C_vs_F] = barrier_potential(Ef, x,
                                                                                    Vc.rows(), Vf.rows(), 0, Vc.rows(),
                                                                                    dhat);
    // compare vertices of C with edges of C
    auto [E_barrier_C_vs_C, g_barrier_C_vs_C, H_barrier_C_vs_C] = barrier_potential(Ec, x,
                                                                                    0, Vc.rows(), 0, Vc.rows(),
                                                                                    dhat);

    // calculating displacement energy and its derivatives manually because TinyAD is annoying when it comes to dynamic
    // size matrices
    double E_disp = (xf - xf_next).squaredNorm() / 2.;
    Eigen::VectorXd g_disp_f = (xf - xf_next);
    Eigen::VectorXd g_disp = vector_segment_assemble(x.size(), xf_start, g_disp_f);
    // don't need to proj E_disp.Hess since its hessian is already SPD
    Eigen::SparseMatrix<double> H_disp(x.size(), x.size());
    std::vector<Eigen::Triplet<double>> list;
    for (long i = xf_start; i < x.size(); i++)
        list.emplace_back(i, i, 1.);
    H_disp.setFromTriplets(list.begin(), list.end());

    double E_fine_inside_coarse = query_mesh_inside(unflatten(xc, 2), Ec, unflatten(xf, 2), Ef)
                                  ? 0 : std::numeric_limits<double>::infinity();

    // energy weights
    double c_cages = 10;
    double c_reinflate = 10000;
    double c_barrier = 10;

    return {c_cages * E_cages
            + c_reinflate * E_disp
            + c_barrier * (E_barrier_F_vs_C + E_barrier_C_vs_F + E_barrier_C_vs_C)
            // following term not differentiable - only takes values {0, inf}
            // only useful for line search, so don't need grad/hess
            + E_fine_inside_coarse,
            c_cages * g_cages
            + c_reinflate * g_disp
            + c_barrier * (g_barrier_F_vs_C + g_barrier_C_vs_F + g_barrier_C_vs_C),
            c_cages * H_cages
            + c_reinflate * H_disp
            + c_barrier * (H_barrier_F_vs_C + H_barrier_C_vs_F + H_barrier_C_vs_C)};
}
