#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <unsupported/Eigen/Splines>
#include <string>
#include <filesystem>
#include "qmoperators/one_electron/AZora/PolyInterpolator.h"

typedef Eigen::Spline<double, 1, 3> Spline1D;
typedef Eigen::SplineFitting<Spline1D> SplineFitting1D;

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
    kappa = kappa * 2.0;
    // The kappa function is half of what is defined in the paper Scalar 
    // Relativistic Effects with Multiwavelets: Implementation and Benchmark
    // it is not used in the code, only the potential is used
}

class RadInterpolater {

    public:
    /**
     * @brief Construct a new Rad Interpolater object
     * @param element The element for which the ZORA potential is to be interpolated
     * @param data_dir The directory containing the ZORA potential data
     * @param mode The mode of interpolation. Either "potential" or "kappa"
     * @param deriv If true, the derivative of the potential is returned
    */
    RadInterpolater(const std::string element, std::string data_dir, const std::string mode, bool deriv=false){
        Eigen::VectorXd rGrid;
        Eigen::VectorXd vZora;
        Eigen::VectorXd kappa;

        this->mode = mode;
        std::string filename = data_dir + '/' + element + ".txt";
        this->deriv = deriv;

        readZoraPotential(filename, rGrid, vZora, kappa);
        if (mode == "kappa") {
            polyZora = std::make_shared<PolyInterpolator>(rGrid, kappa);
        } else if (mode == "potential") {
            polyZora = std::make_shared<PolyInterpolator>(rGrid, vZora);
        } else {
            std::cerr << "Invalid mode. Choose either 'potential' or 'kappa'." << std::endl;
            exit(1);
        }

    }

    double evalf(const double &r) const {
        double y, yp;
        polyZora->evalf(r, y, yp);
        if (deriv) {
            return yp;
        } else {
            return y;
        }
    }

    protected:
    std::shared_ptr<PolyInterpolator> polyZora = nullptr;
    std::string mode;
    bool deriv;

};

// int main() {

//     RadInterpolater spline_V("Ar");

//     Eigen::VectorXd xgrid = Eigen::VectorXd::LinSpaced(100000, 00, 1);
//     Eigen::VectorXd ygrid(xgrid.size());
//     for (int i = 0; i < xgrid.size(); i++) {
//         ygrid(i) = spline_V.evalf(xgrid(i));
//     }

//     // open interpol.txt for writing spline
//     std::ofstream file("interpol.txt");
//     for (int i = 0; i < xgrid.size(); i++) {
//         // write xgrid and ygrid to file using ten digits precision
//         file << std::setprecision(10) << xgrid(i) << " " << ygrid(i) << std::endl;

//     }
//     file.close();

//     return 0;
// }