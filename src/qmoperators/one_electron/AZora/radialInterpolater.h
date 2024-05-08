#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <unsupported/Eigen/Splines>
#include <string>
#include <iostream>
#include <filesystem>

typedef Eigen::Spline<double, 1, 3> Spline1D;
typedef Eigen::SplineFitting<Spline1D> SplineFitting1D;

class RadInterpolater {

    public:
    RadInterpolater(const std::string element){
        Eigen::VectorXd rGrid;
        Eigen::VectorXd vZora;
        Eigen::VectorXd kappa;

        std::filesystem::path p = __FILE__;
        std::filesystem::path parent_dir = p.parent_path();
        std::string filename = parent_dir.string() + element + ".txt";

        readZoraPotential(filename, rGrid, vZora, kappa);
        const auto fitV = SplineFitting1D::Interpolate(kappa.transpose(), 3, rGrid.transpose());
        spline_kappa(fitV);
    }

    double evalf(const double &r){
        return spline_kappa(r).coeff(0);
    }

    protected:
    Spline1D spline_kappa;

}


void readZoraPotential(const std::string path, Eigen::VectorXd &rGrid, Eigen::VectorXd &vZora, Eigen::VectorXd &kappa){
    std::vector<double> r, v, k;
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
}


// int main() {
//     Eigen::VectorXd rGrid;
//     Eigen::VectorXd vZora;
//     Eigen::VectorXd kappa;
//     std::cout << "Reading Zora potential for He" << std::endl;
//     readZoraPotential("He", rGrid, vZora, kappa);
//     std::cout << "Done reading Zora potential for He" << std::endl;
//     std::cout << "size of rGrid: " << rGrid.size() << std::endl;
//     std::cout << "size of vZora: " << vZora.size() << std::endl;


//     const auto fitV = SplineFitting1D::Interpolate(kappa.transpose(), 3, rGrid.transpose());
//     std::cout << "Fitting V" << std::endl;
//     Spline1D spline_V(fitV);
//     std::cout << "Done fitting V" << std::endl;
//     // xgrid is equivalent to np.linspace(0, 10, 10000)
//     Eigen::VectorXd xgrid = Eigen::VectorXd::LinSpaced(100000, 30, 60);
//     Eigen::VectorXd ygrid(xgrid.size());
//     for (int i = 0; i < xgrid.size(); i++) {
//         ygrid(i) = spline_V(xgrid(i)).coeff(0);
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