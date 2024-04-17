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

extern mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;
namespace surface_force {

void plotRandomStuff(mrchem::Molecule &mol, mrchem::OrbitalVector &Phi, double prec){
        auto poisson_pointer = std::make_shared<mrcpp::PoissonOperator>(*mrchem::MRA, prec);

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
    mrcpp::apply(abs_prec, V.real(), *poisson_pointer, rho.real());

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

// Function definition
std::vector<double> surface_forces(mrchem::Molecule &mol, mrchem::OrbitalVector &Phi, double prec) {
    std::cerr << "Surface force calculation" << std::endl;
    plotRandomStuff(mol, Phi, prec);
    std::vector<double> forceValues = {0.0, 0.0, 0.0};
    return forceValues;
}


} // namespace surface_force
