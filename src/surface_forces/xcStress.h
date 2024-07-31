#pragma once

#include "mrchem.h"
#include <Eigen/Core>
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "mrdft/MRDFT.h"
#include <vector>

std::vector<Eigen::Matrix3d> xcLDA(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGrid);
std::vector<Eigen::Matrix3d> xcLDASpin(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGridAlpha, Eigen::MatrixXd &rhoGridBeta);
std::vector<Eigen::Matrix3d> xcGGA(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGrid, Eigen::MatrixXd &nablaRhoGrid);
std::vector<Eigen::Matrix3d> xcGGASpin(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGridAlpha, Eigen::MatrixXd &rhoGridBeta, Eigen::MatrixXd &nablaRhoGridAlpha, Eigen::MatrixXd &nablaRhoGridBeta);