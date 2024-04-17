#include "SurfaceForce.h"

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"

#include "qmoperators/one_electron/NablaOperator.h"
#include <vector>
#include "qmoperators/one_electron/NuclearGradientOperator.h"

#include "chemistry/Molecule.h"
#include "chemistry/Nucleus.h"
#include "chemistry/PhysicalConstants.h"
#include <fstream>
#include <unsupported/Eigen/CXX11/Tensor>

extern mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using namespace Eigen;
using namespace mrchem;
namespace surface_force {

void plotRandomStuff(const mrchem::Molecule &mol, mrchem::OrbitalVector &Phi, double prec){
    auto poisson_op = mrcpp::PoissonOperator(*mrchem::MRA, prec);

    int derivOrder = 1;
    auto mrcd = std::make_shared<mrcpp::BSOperator<3>>(*mrchem::MRA, derivOrder);
    // auto mrcd = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5);
    mrchem::NablaOperator nabla(mrcd);
    nabla.setup(prec);

    // mrcpp::ComplexFunction &rho;
    mrchem::Density rho(false);

    mrchem::density::compute(prec, rho, Phi, DensityType::Total);

    // Adjust precision by system size
    double abs_prec = prec / mrchem::orbital::get_electron_number(Phi);

    mrcpp::ComplexFunction V(false);
    V.alloc(mrchem::NUMBER::Real);
    mrcpp::apply(abs_prec, V.real(), poisson_op, rho.real());

    auto e_field = nabla(V);

    double zmin = -.6;
    double zmax = .6;
    int n = 60;
    auto r = mol.getNuclei()[0].getCoord();

    r[0] = 0.0;
    r[1] = 0.0;
    int Z_k = 1;
    
    for (int i = 0; i < n + 1; i++)
    {
        double z = zmin + i * (zmax - zmin) / (n);
        r[2] = z;
        double c = 0.0001;
        mrchem::NuclearGradientOperator h(Z_k, r, prec, c);
        h.setup(prec);
        auto e = h.trace(Phi).real();
        // std::cerr << "e fielddd " << e[0] << " " << e[1] <<  " " << e[2] << std::endl;
        // std::cerr << "nabla pot " << -o[0].real().evalf(r) << " " << -o[1].real().evalf(r) << " " << -o[2].real().evalf(r) << std::endl;
        std:: cerr << r[2] << " " << e[2] << " " << - e_field[2].real().evalf(r) << " " << V.real().evalf(r) << std::endl;
        h.clear();
    }
    
    // loop over circle in yz plane with radius 0.4 and center at (0,0,0.5)
    double radius = 0.4;

    std::ofstream outfile("toto");
    double pi = 3.14159265359;
    if (outfile.is_open()) {
        for (int i = 0; i < n + 1; i++) {
            r[1] = radius * std::cos(2 * pi * i / n);
            r[2] = radius * std::sin(2 * pi * i / n) + 0.5;
            double c = 0.0001;
            mrchem::NuclearGradientOperator h(Z_k, r, prec, c);
            h.setup(prec);
            auto e = h.trace(Phi).real();
            outfile << 2 * pi * i / n << " " << e[0] <<  " " << e[1] << " " << e[2] << " " << - e_field[0].real().evalf(r) << " " <<  - e_field[1].real().evalf(r) << " " << - e_field[2].real().evalf(r) << " " << std::endl;
            h.clear();
        }
        outfile.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
}


MatrixXd nuclearEfield(const MatrixXd nucPos, const VectorXd nucCharge, const MatrixXd& gridPos) {
    int nGrid = gridPos.rows();
    int nNuc = nucPos.rows();
    MatrixXd Efield = MatrixXd::Zero(nGrid, 3);
    MatrixXd r = MatrixXd::Zero(nGrid, 3);
    double temp;
    for (int i = 0; i < nNuc; ++i) {
        for (int j = 0; j < nGrid; j++)
        {
            r.row(j) = nucPos.row(i) - gridPos.row(j);
            temp = r.row(j).norm();
            temp = temp * temp * temp;
            Efield.row(j) += nucCharge(i) * r.row(j) / temp;
        }
    }
    return Efield;
}

mrcpp::ComplexFunction calcPotential(Density rho, mrcpp::PoissonOperator poisson, double prec) {
    mrcpp::ComplexFunction V(false);
    V.alloc(mrchem::NUMBER::Real);
    mrcpp::apply(prec, V.real(), poisson, rho.real());
    return V;
}

MatrixXd electronicEfield(mrcpp::ComplexFunction pot, NablaOperator nabla, const MatrixXd& gridPos) {
    int nGrid = gridPos.rows();
    MatrixXd Efield = MatrixXd::Zero(nGrid, 3);
    auto e_field = nabla(pot);
    for (int i = 0; i < nGrid; i++)
    {
        std::array<double, 3> pos = {gridPos(i, 0), gridPos(i, 1), gridPos(i, 2)};
        Efield(i, 0) = -e_field[0].real().evalf(pos);
        Efield(i, 1) = -e_field[1].real().evalf(pos);
        Efield(i, 2) = -e_field[2].real().evalf(pos);
    }
    return Efield;
}

Eigen::Tensor<double, 3> maxwellStress(const Molecule mol, mrcpp::ComplexFunction pot, NablaOperator nabla, const MatrixXd gridPos){
    int nGrid = gridPos.rows();
    int nNuc = mol.getNNuclei();

    Eigen::MatrixXd nucPos(nNuc, 3);
    Eigen::VectorXd nucCharge(nNuc);
    for (int i = 0; i < nNuc; i++)
    {
        std::array<double, 3> coord = mol.getNuclei()[i].getCoord();
        nucPos(i, 0) = coord[0];
        nucPos(i, 1) = coord[1];
        nucPos(i, 2) = coord[2];
        nucCharge(i) = mol.getNuclei()[i].getCharge();
    }

    MatrixXd Efield = electronicEfield(pot, nabla, gridPos) + nuclearEfield(nucPos, nucCharge, gridPos);

    Eigen::Tensor<double, 3> stress(nGrid, 3, 3);
    for (int i = 0; i < nNuc; i++) {
        for (int i1 = 0; i1 < 3; i1++) {
            for (int i2 = 0; i2 < 3; i2++) {
                stress(i, i1, i2) = Efield(i, i1) * Efield(i, i2);
            }
        }
        for (int i1 = 0; i1 < 3; i1++){
            stress(i, i1, i1) -= 0.5 * Efield(i, 0) * Efield(i, 0) + Efield(i, 1) * Efield(i, 1) + Efield(i, 2) * Efield(i, 2);
        }
        for (int i1 = 0; i1 < 3; i1++){
            for (int i2 = 0; i2 < 3; i2++){
                stress(i, i1, i2) *= 1.0 / (4 * M_PI);
            }
        }   
    }
    return stress;
}


// Function definition
std::vector<double> surface_forces(mrchem::Molecule &mol, mrchem::OrbitalVector &Phi, double prec) {
    std::cerr << "Surface force calculation" << std::endl;
    plotRandomStuff(mol, Phi, prec);
    std::vector<double> forceValues = {0.0, 0.0, 0.0};
    return forceValues;
}


} // namespace surface_force
