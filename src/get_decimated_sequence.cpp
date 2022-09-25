#include "get_decimated_sequence.h"

bool get_decimated_sequence(Eigen::MatrixXd &V, Eigen::MatrixXi &F, int num_meshes, int coarsest_faces,
                            std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> &out)
{
    for (int i = 0; i < num_meshes; i++)
    {
        double t = (double) (i + 1) / num_meshes;
        int num_faces = F.rows() * (1 - t) + coarsest_faces * t;
        Eigen::MatrixXd U;
        Eigen::MatrixXi G;
        Eigen::VectorXi J;
        if (not igl::decimate(V, F, num_faces, U, G, J))
            return false;
        out.emplace_back(U, G);
    }
    return true;
}