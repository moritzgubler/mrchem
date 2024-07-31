#pragma once

#include "mrchem.h"
#include <Eigen/Core>
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "mrdft/MRDFT.h"
#include <vector>

std::vector<Eigen::Matrix3d> xcLDAStress(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGrid);
std::vector<Eigen::Matrix3d> xcLDASpinStress(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGridAlpha, Eigen::MatrixXd &rhoGridBeta);
std::vector<Eigen::Matrix3d> xcGGAStress(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGrid, Eigen::MatrixXd &nablaRhoGrid);
std::vector<Eigen::Matrix3d> xcGGASpinStress(std::unique_ptr<mrdft::MRDFT> &mrdft_p, Eigen::MatrixXd &rhoGridAlpha, Eigen::MatrixXd &rhoGridBeta, Eigen::MatrixXd &nablaRhoGridAlpha, Eigen::MatrixXd &nablaRhoGridBeta);