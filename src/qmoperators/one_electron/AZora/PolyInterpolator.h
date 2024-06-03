#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <filesystem>
#include <iostream>

/**
 * @brief Read ZORA potential from file. Check if file exists and abort if it does not.
 * @param path Path to the file containing the ZORA potential
 * @param rGrid Vector containing the radial grid
 * @param vZora Vector containing the ZORA potential
 * @param kappa Vector containing the kappa parameter
*/
void readZoraPotential(const std::string path, Eigen::VectorXd &rGrid, Eigen::VectorXd &vZora, Eigen::VectorXd &kappa){
    std::vector<double> r, v, k;
    bool file_exists = std::filesystem::exists(path);
    if (!file_exists) {
        std::cerr << "File " << path << " does not exist." << std::endl;
        std::cout << "File " << path << " does not exist." << std::endl;
        exit(1);
    }
    std::ifstream file(path);
    std::string line;
    double r_, v_, k_;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        iss >> k_ >> r_ >> v_;
        r.push_back(r_);
        v.push_back(v_);
        k.push_back(k_);
    }
    file.close();
    rGrid = Eigen::Map<Eigen::VectorXd>(r.data(), r.size());
    vZora = Eigen::Map<Eigen::VectorXd>(v.data(), v.size());
    kappa = Eigen::Map<Eigen::VectorXd>(k.data(), k.size());
    // The kappa function is half of what is defined in the paper Scalar 
    // Relativistic Effects with Multiwavelets: Implementation and Benchmark
    // it is not used in the code, only the potential is used
}

double polynomialInterpolate5(Eigen::VectorXd &x_in, Eigen::VectorXd &y_in, double x){
    double xm2 = x_in(0);
    double xm1 = x_in(1);
    double x00 = x_in(2);
    double xp1 = x_in(3);
    double xp2 = x_in(4);

    double ym2 = y_in(0);
    double ym1 = y_in(1);
    double y00 = y_in(2);
    double yp1 = y_in(3);
    double yp2 = y_in(4);
    return ym2 + (x - xm2)*((ym1 - ym2)/(xm1 - xm2) +
        (x - xm1)*(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2) +
        (x - x00)*((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
        ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1) +
        ((x - xp1)*(-((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
        ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/
        (-xm2 + xp1)) + (-(((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/
        (-xm1 + xp1)) + ((-y00 + yp1)/(x00 - xp1) + (yp1 - yp2)/(xp1 - xp2))/(-x00 + xp2))/
        (-xm1 + xp2)))/(-xm2 + xp2))));
}

double polynomialInterpolate5_deriv(Eigen::VectorXd &x_in, Eigen::VectorXd &y_in, double x) {
    double xm2 = x_in(0);
    double xm1 = x_in(1);
    double x00 = x_in(2);
    double xp1 = x_in(3);
    double xp2 = x_in(4);

    double ym2 = y_in(0);
    double ym1 = y_in(1);
    double y00 = y_in(2);
    double yp1 = y_in(3);
    double yp2 = y_in(4);

    return (ym1 - ym2)/(xm1 - xm2) + (x - xm1)*(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2) +
             (x - x00)*((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
                 ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1) +
              ((x - xp1)*(-((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
                        ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1)) +
                   (-(((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1)) +
                      ((-y00 + yp1)/(x00 - xp1) + (yp1 - yp2)/(xp1 - xp2))/(-x00 + xp2))/(-xm1 + xp2)))/(-xm2 + xp2))) +
        (x - xm2)*(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2) +
           (x - x00)*((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
                 ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1) +
              ((x - xp1)*(-((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
                        ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1)) +
                   (-(((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1)) +
                      ((-y00 + yp1)/(x00 - xp1) + (yp1 - yp2)/(xp1 - xp2))/(-x00 + xp2))/(-xm1 + xp2)))/(-xm2 + xp2)) +
           (x - xm1)*((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
                 ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1) +
              ((x - x00)*(-((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
                        ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1)) +
                   (-(((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1)) +
                      ((-y00 + yp1)/(x00 - xp1) + (yp1 - yp2)/(xp1 - xp2))/(-x00 + xp2))/(-xm1 + xp2)))/(-xm2 + xp2) +
              ((x - xp1)*(-((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
                        ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1)) +
                   (-(((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1)) +
                      ((-y00 + yp1)/(x00 - xp1) + (yp1 - yp2)/(xp1 - xp2))/(-x00 + xp2))/(-xm1 + xp2)))/(-xm2 + xp2)));
}

class PolyInterpolator {
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    double xmin;
    double xmax;
    public:
    PolyInterpolator(Eigen::VectorXd &x_in, Eigen::VectorXd &y_in){
        x = x_in;
        y = y_in;
        xmin = x_in(0);
        xmax = x_in(x_in.size() - 1);
    }
    
    void evalf(const double &x, double &y, double &yp) const {
        if (x < xmin || x > xmax) {
            if (x < xmin) {
                y = this->y(0);
                yp = 0;
            } else {
                y = this->y(this->y.size() - 1);
                yp = 0;
            }
        }

        // Find the interval in which x lies using binary search
        int i = 0;
        int j = this->x.size() - 1;
        while (j - i > 1) {
            int k = (i + j)/2;
            if (x < this->x(k)) {
                j = k;
            } else {
                i = k;
            }
        }

        Eigen::VectorXd x_in(5), y_in(5);

        if (i == 0) i = 2;
        if (i == 1) i = 2;
        if (i == this->x.size() - 1) i = this->x.size() - 3;
        if (i == this->x.size() - 2) i = this->x.size() - 3;
        for (int k = 0; k < 5; k++) {
            x_in(k) = this->x(i - 2 + k);
            y_in(k) = this->y(i - 2 + k);
        }

        y = polynomialInterpolate5(x_in, y_in, x);
        yp = polynomialInterpolate5_deriv(x_in, y_in, x);

    }
};