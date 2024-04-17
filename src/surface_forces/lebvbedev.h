#ifndef LEBVEDEV_H
#define LEBVEDEV_H

#include "lebvbedev.h"
#include <Eigen/Dense>
#include <string>

class LebedevIntegrator {
public:
    Eigen::MatrixXd points;
    Eigen::VectorXd weights;
    Eigen::MatrixXd normals;
    int n;


    LebedevIntegrator(const std::string& filename, double radius, const Eigen::Vector3d& center);
    Eigen::MatrixXd getPoints() const;
    Eigen::VectorXd getWeights() const;
    Eigen::MatrixXd getNormals() const;
};

#endif // LEBVEDEV_H