#include "helpers.h"

Eigen::VectorXd flatten(const Eigen::MatrixXd &V) {
    Eigen::Index rows = V.rows(), cols = V.cols();
    Eigen::VectorXd out = Eigen::VectorXd::Zero(rows * cols);
    for (int r = 0; r < rows; r++) {
        out.segment(r * cols, cols) = V.row(r);
    }
    return out;
}

Eigen::MatrixXd unflatten(const Eigen::VectorXd &v, int dim) {
    assert(dim > 0);
    assert(v.rows() % dim == 0);
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(v.rows() / dim, dim);
    for (int i = 0; i < v.rows(); i++) {
        out(i / dim, i % dim) = v(i);
    }
    return out;
}

Eigen::VectorXd vector_segment_assemble(Eigen::Index size, Eigen::Index start, const Eigen::VectorXd &segment) {
    Eigen::VectorXd out = Eigen::VectorXd::Zero(size);
    out.segment(start, segment.rows()) = segment;
    return out;
}

Eigen::SparseMatrix<double> sparse_block_assemble(Eigen::Index rows, Eigen::Index cols,
                                                  Eigen::Index startRow, Eigen::Index startCol,
                                                  const Eigen::SparseMatrix<double> &block) {
    Eigen::SparseMatrix<double> out(rows, cols);
    std::vector<Eigen::Triplet<double>> list;
    for (int i = 0; i < block.outerSize(); i++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(block, i); it; ++it) {
            list.emplace_back(it.row() + startRow, it.col() + startCol, it.value());
        }
    }
    out.setFromTriplets(list.begin(), list.end());
    return out;
}

int query_winding_number(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E,
                         const Eigen::RowVector2d &p) {
    // TODO: not sure if robust implementation
    Eigen::Matrix2d A;
    A.col(0) = Eigen::Vector2d(1, 0);
    int winding = 0;
    for (int e = 0; e < E.rows(); e++) {
        int vc0 = E(e, 0), vc1 = E(e, 1);
        A.col(1) = (V.row(vc1) - V.row(vc0)).transpose();
        Eigen::Vector2d c = p - V.row(vc0);
        if (A.determinant() != 0) {
            Eigen::Vector2d st = A.inverse() * c;
            if (st(0) < 0. and 0. < st(1) and st(1) < 1.)
                winding++;
        }
    }
    return winding % 2;
}

bool query_point_inside(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::Vector2d &p) {
    return query_winding_number(V, E, p) % 2 == 1;
}

bool query_mesh_inside(const Eigen::MatrixXd &V_outside, const Eigen::MatrixXi &E_outside,
                       const Eigen::MatrixXd &V_inside, const Eigen::MatrixXi &E_inside) {
    for (int v = 0; v < V_inside.rows(); v++)
        if (!query_point_inside(V_outside, E_outside, V_inside.row(v)))
            return false;
    return !query_meshes_intersect(flatten(V_outside), E_outside, flatten(V_inside), E_inside);
}

bool query_meshes_intersect(const Eigen::VectorXd &x0, const Eigen::MatrixXi &E0,
                            const Eigen::VectorXd &x1, const Eigen::MatrixXi &E1) {
    for (int e0 = 0; e0 < E0.rows(); e0++) {
        Eigen::Vector2d a0 = x0.segment<2>(E0(e0, 0) * 2);
        Eigen::Vector2d d0 = x0.segment<2>(E0(e0, 1) * 2) - a0;
        for (int e1 = 0; e1 < E1.rows(); e1++) {
            Eigen::Vector2d a1 = x1.segment<2>(E1(e1, 0) * 2);
            Eigen::Vector2d d1 = x1.segment<2>(E1(e1, 1) * 2) - a1;
            Eigen::Matrix2d A;
            A.col(0) = d0;
            A.col(1) = -d1;
            if (A.determinant() == 0) {
                continue;
                // TODO: check if edges overlap
            }
            Eigen::Vector2d st = A.inverse() * (a1 - a0);
            double s = st(0), t = st(1);
            if (0 <= s and s <= 1 and 0 <= t and t <= 1)
                return true;
        }
    }
    return false;
}
