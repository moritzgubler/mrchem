#include "surface_forces/xcStress.h"

#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "mrdft/MRDFT.h"

using namespace Eigen;

std::vector<Eigen::Matrix3d> xcLDA(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGrid){
    int nGrid = rhoGrid.rows();
    std::vector<Eigen::Matrix3d> out(nGrid);
    Eigen::MatrixXd xcOUT =  mrdft_p->functional().evaluate_transposed(rhoGrid);
    for (int i = 0; i < nGrid; i++) {
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
    std::vector<Matrix3d> out = std::vector<Eigen::Matrix3d>(nGrid);
    inp.col(0) = rhoGridAlpha.col(0);
    inp.col(1) = rhoGridBeta.col(0);
    Eigen::MatrixXd xc = mrdft_p->functional().evaluate_transposed(inp);
    for (int i = 0; i < rhoGridAlpha.rows(); i++) {
        out[i] = Matrix3d::Zero();
        for (int j = 0; j < 3; j++) {
            out[i](j, j) = xc(i, 0) - xc(i, 1) * rhoGridAlpha(i) - xc(i, 2) * rhoGridBeta(i);
        }
    }
    return out;
}

std::vector<Matrix3d> xcGGA(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGrid, Eigen::MatrixXd &nablaRhoGrid){
    int nGrid = rhoGrid.rows();
    Eigen::MatrixXd inp(rhoGrid.rows(), 4);
    inp.col(0) = rhoGrid.col(0);
    inp.col(1) = nablaRhoGrid.col(0);
    inp.col(2) = nablaRhoGrid.col(1);
    inp.col(3) = nablaRhoGrid.col(2);
    Eigen::MatrixXd xcOUT =  mrdft_p->functional().evaluate_transposed(inp);
    std::vector<Matrix3d> out(nGrid);
    for (int i = 0; i < rhoGrid.rows(); i++) {
        out[i] = Matrix3d::Zero();

        for (int j = 0; j < 3; j++) {
            out[i](j, j) = xcOUT(i, 0) - xcOUT(i, 1) * rhoGrid(i);
        }
        for (int j1 = 0; j1 < 3; j1++) {
            for (int j2 = 0; j2 < 3; j2++) {
                out[i](j1, j2) = out[i](j1, j2) - xcOUT(i, 2 + j1) * nablaRhoGrid(i, j2) * rhoGrid(i);
            }
        }
    }
    return out;
}

std::vector<Matrix3d> xcGGASpin(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGridAlpha, Eigen::MatrixXd &rhoGridBeta, Eigen::MatrixXd &nablaRhoGridAlpha, Eigen::MatrixXd &nablaRhoGridBeta){
    int nGrid = rhoGridAlpha.rows();
    Eigen::MatrixXd inp(rhoGridAlpha.rows(), 8);
    std::vector<Matrix3d> out = std::vector<Eigen::Matrix3d>(nGrid);
    inp.col(0) = rhoGridAlpha.col(0);
    inp.col(1) = rhoGridBeta.col(0);
    inp.col(2) = nablaRhoGridAlpha.col(0);
    inp.col(3) = nablaRhoGridAlpha.col(1);
    inp.col(4) = nablaRhoGridAlpha.col(2);
    inp.col(5) = nablaRhoGridBeta.col(0);
    inp.col(6) = nablaRhoGridBeta.col(1);
    inp.col(7) = nablaRhoGridBeta.col(2);
    Eigen::MatrixXd xc = mrdft_p->functional().evaluate_transposed(inp);
    for (int i = 0; i < rhoGridAlpha.rows(); i++) {
        out[i] = Matrix3d::Zero();
        for (int j = 0; j < 3; j++) {
            out[i](j, j) = xc(i, 0) - xc(i, 1) * rhoGridAlpha(i) - xc(i, 2) * rhoGridBeta(i);
        }
        for (int j1 = 0; j1 < 3; j1++) {
            for (int j2 = 0; j2 < 3; j2++) {
                out[i](j1, j2) = out[i](j1, j2) 
                    - xc(i, 3 + j1) * nablaRhoGridAlpha(i, j2) * rhoGridAlpha(i) 
                    - xc(i, 6 + j1) * nablaRhoGridBeta(i, j2) * rhoGridBeta(i);
            }
        }
    }
    return out;
}