#pragma once

#include "mrchem.h"
#include <Eigen/Core>
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"

Eigen::MatrixXd createXCInput(mrchem::Density &rho, mrchem::OrbitalVector &nablaRho, Eigen::MatrixXd &gridPos);