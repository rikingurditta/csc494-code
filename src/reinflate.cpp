#include "reinflate.h"

void reinflate(const Eigen::MatrixXd &Vc, const Eigen::MatrixXi &Ec,
               const Eigen::MatrixXd &Vf, const Eigen::MatrixXi &Ef,
               const std::vector<Eigen::MatrixXd> &Vf_flow,
               Eigen::MatrixXd &Vc_inflated)
{
    Vc_inflated = Vc;
    Eigen::MatrixXd Vc_curr = Vc;
    for (int i = 0; i < Vf_flow.size() - 1; i++)
    {
        Eigen::MatrixXd Vf_curr = Vf_flow[Vf_flow.size() - 1 - i];
        Eigen::MatrixXd Vf_next = Vf_flow[Vf_flow.size() - 1 - i - 1];
        // inflate current coarse mesh by one flow step
        reinflate_simulate_timestep(Vc_curr, Ec, Vf_curr, Vf_next, Ef, Vc_inflated);
        Vc_curr = Vc_inflated;
    }
}

void reinflate_simulate_timestep(const Eigen::MatrixXd &Vc, const Eigen::MatrixXi &Ec,
                                 const Eigen::MatrixXd &Vf_curr, const Eigen::MatrixXd &Vf_next, const Eigen::MatrixXi &Ef,
                                 Eigen::MatrixXd &Vc_inflated)
{
}