#include <igl/opengl/glfw/Viewer.h>

#include "get_decimated_sequence.h"
#include "flow.h"
#include "igl/readCSV.h"
#include "reinflate.h"

int main(int argc, char *argv[])
{
    Eigen::MatrixXd Vf, Vc;
    Eigen::MatrixXi Ef, Ec;
    igl::readCSV("../data/2d_torus_V.csv", Vf);
    igl::readCSV("../data/2d_torus_E.csv", Ef);
    igl::readCSV("../data/2d_torus_dec1_V.csv", Vc);
    igl::readCSV("../data/2d_torus_dec1_E.csv", Ec);

    std::vector<Eigen::MatrixXd> Vf_flow;
    flow(Vc, Ec, Vf, Ef, Vf_flow);
    Eigen::MatrixXd Vc_inflated;
//    reinflate(Vc, Ec, Vf, Ef, Vf_flow, Vc_inflated);
}
