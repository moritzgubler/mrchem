#pragma once

#include "mrchem.h"
#include <Eigen/Core>
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "mrdft/MRDFT.h"

Eigen::MatrixXd xcLDA(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGrid);
Eigen::MatrixXd xcLDASpin(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGridAlpha, Eigen::MatrixXd &rhoGridBeta);
Eigen::MatrixXd createXCInput(mrchem::Density &rho, mrchem::OrbitalVector &nablaRho, Eigen::MatrixXd &gridPos);