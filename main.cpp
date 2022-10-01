#include <igl/opengl/glfw/Viewer.h>

#include "get_decimated_sequence.h"
#include "flow.h"
#include "reinflate.h"
#include "igl/readCSV.h"

int main(int argc, char *argv[])
{
    Eigen::MatrixXd V;
    igl::readCSV("../data/2d_mesh_V.csv", V);
    Eigen::MatrixXi E;
    igl::readCSV("../data/2d_mesh_E.csv", E);

}
