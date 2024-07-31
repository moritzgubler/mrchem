#include "surface_forces/xcHelper.h"

#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "mrdft/MRDFT.h"

using namespace Eigen;

std::vector<Eigen::Matrix3d> xcLDA(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGrid){
    int nGrid = rhoGrid.rows();
    std::vector<Eigen::Matrix3d> out(nGrid);
    // Eigen::MatrixXd out = Eigen::MatrixXd::Zero(rhoGrid.rows(), 5);
    Eigen::MatrixXd xcOUT =  mrdft_p->functional().evaluate_transposed(rhoGrid);
    // out.col(0) = xcOUT.col(0);
    for (int i = 0; i < nGrid; i++) {
        // out(i, 1) = xcOUT(i, 1) * rhoGrid(i);
        out[i] = Matrix3d::Zero();
        for (int j = 0; j < 3; j++) {
            out[i](j, j) = xcOUT(i, 0) - xcOUT(i, 1) * rhoGrid(i);
        }
    }

    return out;
}

std::vector<Matrix3d> xcLDASpin(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGridAlpha, Eigen::MatrixXd &rhoGridBeta){
    int nGrid = rhoGridAlpha.rows();
    Eigen::MatrixXd inp(rhoGridAlpha.rows(), 2);
    // Eigen::MatrixXd out = Eigen::MatrixXd::Zero(rhoGridAlpha.rows(), 5);
    std::vector<Matrix3d> out = std::vector<Eigen::Matrix3d>(nGrid);
    // Eigen::MatrixXd outAlpha = Eigen::MatrixXd::Zero(rhoGridAlpha.rows(), 5);
    // Eigen::MatrixXd outBeta = Eigen::MatrixXd::Zero(rhoGridAlpha.rows(), 5);
    inp.col(0) = rhoGridAlpha.col(0);
    inp.col(1) = rhoGridBeta.col(0);
    Eigen::MatrixXd xc = mrdft_p->functional().evaluate_transposed(inp);
    // out.col(0) = xc.col(0);
    for (int i = 0; i < rhoGridAlpha.rows(); i++) {
        // out(i, 1) = xc(i, 1) * rhoGridAlpha(i) + xc(i, 2) * rhoGridBeta(i);
        out[i] = Matrix3d::Zero();
        for (int j = 0; j < 3; j++) {
            out[i](j, j) = xc(i, 0) - xc(i, 1) * rhoGridAlpha(i) - xc(i, 2) * rhoGridBeta(i);
        }
    }
    return out;
}

// Eigen::MatrixXd xcGGA(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGrid, Eigen::MatrixXd &nablaRhoGrid){
//     Eigen::MatrixXd out = Eigen::MatrixXd::Zero(rhoGrid.rows(), 5);
//     Eigen::MatrixXd inp(rhoGrid.rows(), 4);
//     inp.col(0) = rhoGrid.col(0);
//     inp.col(1) = nablaRhoGrid.col(0);
//     inp.col(2) = nablaRhoGrid.col(1);
//     inp.col(3) = nablaRhoGrid.col(2);
//     Eigen::MatrixXd xcOUT =  mrdft_p->functional().evaluate_transposed(inp);
//     out.col(0) = xcOUT.col(0);
//     for (int i = 0; i < rhoGrid.rows(); i++) {
//         out(i, 1) = xcOUT(i, 1) * rhoGrid(i);
//         out(i, 2) = xcOUT(i, 2) * nablaRhoGrid(i, 0);
//         out(i, 3) = xcOUT(i, 3) * nablaRhoGrid(i, 1);
//         out(i, 4) = xcOUT(i, 4) * nablaRhoGrid(i, 2);
//     }

//     return out;
// }

