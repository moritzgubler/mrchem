#include "surface_forces/xcHelper.h"

#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "mrdft/MRDFT.h"

Eigen::MatrixXd xcLDA(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGrid){
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(rhoGrid.rows(), 5);
    Eigen::MatrixXd xcOUT =  mrdft_p->functional().evaluate_transposed(rhoGrid);
    out.col(0) = xcOUT.col(0);
    for (int i = 0; i < rhoGrid.rows(); i++) {
        out(i, 1) = xcOUT(i, 1) * rhoGrid(i);
    }

    return out;
}

Eigen::MatrixXd xcLDASpin(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGridAlpha, Eigen::MatrixXd &rhoGridBeta){
    Eigen::MatrixXd inp(rhoGridAlpha.rows(), 2);
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(rhoGridAlpha.rows(), 5);
    Eigen::MatrixXd outAlpha = Eigen::MatrixXd::Zero(rhoGridAlpha.rows(), 5);
    Eigen::MatrixXd outBeta = Eigen::MatrixXd::Zero(rhoGridAlpha.rows(), 5);
    inp.col(0) = rhoGridAlpha.col(0);
    inp.col(1) = rhoGridBeta.col(0);
    Eigen::MatrixXd xc = mrdft_p->functional().evaluate_transposed(inp);
    out.col(0) = xc.col(0);
    for (int i = 0; i < rhoGridAlpha.rows(); i++) {
        out(i, 1) = xc(i, 1) * rhoGridAlpha(i) + xc(i, 2) * rhoGridBeta(i);
    }
    
    return out;
}

