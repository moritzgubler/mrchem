#include "surface_forces/xcHelper.h"

#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "mrdft/MRDFT.h"


Eigen::MatrixXd xcLDA(std::shared_ptr<mrdft::MRDFT> mrdft_p, Eigen::MatrixXd &rhoGrid){
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(5, rhoGrid.cols());
    out(Eigen::seq(0, 1), Eigen::all) =  mrdft_p->functional().evaluate_transposed(rhoGrid);
    return out;
}

Eigen::MatrixXd xcLDASpin(std::shared_ptr<mrdft::MRDFT> mrdft_p, Eigen::MatrixXd &rhoGridAlpha, Eigen::MatrixXd &rhoGridBeta){
    Eigen::MatrixXd inp(2, rhoGridAlpha.cols());
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(4, rhoGridAlpha.cols());
    inp.row(0) = rhoGridAlpha;
    inp.row(1) = rhoGridBeta;
    Eigen::MatrixXd xc = mrdft_p->functional().evaluate_transposed(inp);
    out.row(0) = xc.row(0);
    out.row(1) = xc.row(1) + xc.row(2);
    return out;
}

Eigen::MatrixXd createXCInput(mrchem::Density &rho, mrchem::OrbitalVector &nablaRho, Eigen::MatrixXd &gridPos){

}