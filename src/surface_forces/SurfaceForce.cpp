#include "surface_forces/SurfaceForce.h"

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

#include "mrdft/Factory.h"
#include "mrdft/MRDFT.h"
#include "qmoperators/two_electron/XCOperator.h"
#include "mrdft/Functional.h"

#include "qmoperators/one_electron/HessianOperator.h"
#include "tensor/RankOneOperator.h"

#include "surface_forces/lebvedev.h"
#include <string>
#include <iostream>
#include <filesystem>

extern mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using namespace Eigen;
using namespace mrchem;
using nlohmann::json;

namespace surface_force {

MatrixXd nuclearEfield(const MatrixXd &nucPos, const VectorXd &nucCharge, const MatrixXd gridPos) {
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

mrcpp::ComplexFunction calcPotential(Density &rho, mrcpp::PoissonOperator &poisson, double prec) {
    mrcpp::ComplexFunction V(false);
    V.alloc(mrchem::NUMBER::Real);
    mrcpp::apply(prec, V.real(), poisson, rho.real());
    return V;
}

MatrixXd electronicEfield(mrcpp::ComplexFunction &pot, NablaOperator &nabla, const MatrixXd &gridPos) {
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

std::vector<Eigen::Matrix3d> maxwellStress(const Molecule &mol, mrcpp::ComplexFunction &pot, NablaOperator &nabla, const MatrixXd &gridPos){
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

    std::vector<Eigen::Matrix3d> stress(nGrid);
    for (int i = 0; i < nGrid; i++) {
        for (int i1 = 0; i1 < 3; i1++) {
            for (int i2 = 0; i2 < 3; i2++) {
                stress[i](i1, i2) = Efield(i, i1) * Efield(i, i2);
            }
        }
        for (int i1 = 0; i1 < 3; i1++){
            stress[i](i1, i1) = stress[i](i1, i1) - 0.5 * (Efield(i, 0) * Efield(i, 0) + Efield(i, 1) * Efield(i, 1) + Efield(i, 2) * Efield(i, 2));
        }
        for (int i1 = 0; i1 < 3; i1++){
            for (int i2 = 0; i2 < 3; i2++){
                stress[i](i1, i2) *= 1.0 / (4 * M_PI);
            }
        }
    }
    return stress;
}

std::vector<Matrix3d> xcStress(const Molecule &mol, const Density &rho, const MatrixXd &gridPos, const json &json_fock, double prec){
    int nGrid = gridPos.rows();

    int order = 0;


    std::vector<Matrix3d> stress(nGrid);

    Eigen::MatrixXd rhoGrid(nGrid, 1);
    std::array<double, 3> pos;
    for (int i = 0; i < nGrid; i++)
    {
        pos[0] = gridPos(i, 0);
        pos[1] = gridPos(i, 1);
        pos[2] = gridPos(i, 2);
        rhoGrid(i, 0) = rho.real().evalf(pos);

    }
    // std::cerr << "rhoGrid" << std::endl;

    bool shared_memory = json_fock["xc_operator"]["shared_memory"];
    auto json_xcfunc = json_fock["xc_operator"]["xc_functional"];
    auto xc_spin = json_xcfunc["spin"];
    auto xc_cutoff = json_xcfunc["cutoff"];
    auto xc_funcs = json_xcfunc["functionals"];
    auto xc_order = order + 1;

    auto Phi_p = mol.getOrbitals_p();

    mrdft::Factory xc_factory(*MRA);
    xc_factory.setSpin(xc_spin);
    xc_factory.setOrder(xc_order);
    xc_factory.setDensityCutoff(xc_cutoff);
    for (const auto &f : xc_funcs) {
        auto name = f["name"];
        auto coef = f["coef"];
        xc_factory.setFunctional(name, coef);
    }

    std::unique_ptr<mrdft::MRDFT> mrdft_p = xc_factory.build();
    auto XC_p = std::make_shared<XCOperator>(mrdft_p, Phi_p, shared_memory);

    XC_p->potential->setup(prec);
    VectorXd xcGrid(nGrid);
    VectorXd vxcGrid(nGrid);
    // open file:
    // std::ofstream outfile("toto_xc");
    for (int i = 0; i < nGrid; i++){
        pos[0] = gridPos(i, 0);
        pos[1] = gridPos(i, 1);
        pos[2] = gridPos(i, 2);
        xcGrid(i) = std::get<1>(XC_p->potential->xc_vec[0])->evalf(pos);
        vxcGrid(i) = std::get<1>(XC_p->potential->xc_vec[1])->evalf(pos);
        // std::cerr << pos[2] << " " << rhoGrid(i) << " " << xcGrid(i) << " " << vxcGrid(i) << std::endl;
        // outfile << pos[2] << " " << rhoGrid(i) << " " << xcGrid(i) << " " << vxcGrid(i) << std::endl;
        for (int i1 = 0; i1 < 3; i1++) {
            for (int i2 = 0; i2 < 3; i2++) {
                stress[i](i1, i2) = 0.0;
            }
        }
        for (int i1 = 0; i1 < 3; i1++)
        {
            stress[i](i1, i1) = xcGrid(i) - vxcGrid(i) * rhoGrid(i);
        }
        
    }
    // outfile.close();
    
    XC_p->potential->clear();

    return stress;
}

std::vector<Matrix3d> kineticStress(const mrchem::Molecule &mol, mrchem::OrbitalVector &Phi, double prec, const MatrixXd &gridPos){

    int nGrid = gridPos.rows();
    int nOrbs = Phi.size();

    int derivOrder1 = 1;
    auto D1 = std::make_shared<mrcpp::BSOperator<3>>(*mrchem::MRA, derivOrder1);
    int derivOrder2 = 2;
    auto D2 = std::make_shared<mrcpp::BSOperator<3>>(*mrchem::MRA, derivOrder2);
    mrchem::HessianOperator hess(D1, D2, prec);
    hess.setup(prec);
    mrchem::NablaOperator nabla(D1);
    nabla.setup(prec);

    double orbVal;
    std::vector<Matrix3d> stress(nGrid);
    for (int i = 0; i < nGrid; i++) {
        stress[i] = Matrix3d::Zero();
    }
    std::array<double, 3> pos;
    double n1, n2, n3;
    for (int iOrb = 0; iOrb < Phi.size(); iOrb++) {
    
        mrcpp::ComplexFunction orb = Phi[iOrb];
        auto hessPhi = hess(orb);
        auto nablaPhi = nabla(orb);

        for (int i = 0; i < nGrid; i++) {
            pos[0] = gridPos(i, 0);
            pos[1] = gridPos(i, 1);
            pos[2] = gridPos(i, 2);
            orbVal = Phi[iOrb].real().evalf(pos);
            n1 = nablaPhi[0].real().evalf(pos);
            n2 = nablaPhi[1].real().evalf(pos);
            n3 = nablaPhi[2].real().evalf(pos);
            stress[i](0, 0) += orbVal * hessPhi[0].real().evalf(pos) - n1 * n1;
            stress[i](1, 1) += orbVal * hessPhi[1].real().evalf(pos) - n2 * n2;
            stress[i](2, 2) += orbVal * hessPhi[2].real().evalf(pos) - n3 * n3;
            stress[i](0, 1) += orbVal * hessPhi[3].real().evalf(pos) - n1 * n2;
            stress[i](0, 2) += orbVal * hessPhi[4].real().evalf(pos) - n1 * n3;
            stress[i](1, 2) += orbVal * hessPhi[5].real().evalf(pos) - n2 * n3;
            stress[i](1, 0) = stress[i](0, 1);
            stress[i](2, 0) = stress[i](0, 2);
            stress[i](2, 1) = stress[i](1, 2);
        }
    }
    nabla.clear();
    hess.clear();
    return stress;
}

/**
 * Calculates the forces using surface integrals for a given molecule and orbital vector.
 *
 * @param mol The molecule for which to calculate the forces.
 * @param Phi The orbital vector obtained from the SCF calculation.
 * @param prec The precision value used in the calculation.
 * @param json_fock The JSON object containing the Fock matrix settings.
 * @return The matrix of forces, shape (nAtoms, 3).
 */
Eigen::MatrixXd surface_forces(mrchem::Molecule &mol, mrchem::OrbitalVector &Phi, double prec, const json &json_fock) {

    // setup density
    mrchem::Density rho(false);
    mrchem::density::compute(prec, rho, Phi, DensityType::Total);

    // setup operators and potentials
    auto poisson_op = mrcpp::PoissonOperator(*mrchem::MRA, prec);
    int derivOrder = 1;
    auto mrcd = std::make_shared<mrcpp::BSOperator<3>>(*mrchem::MRA, derivOrder);
    mrchem::NablaOperator nabla(mrcd);
    nabla.setup(prec);
    double abs_prec = prec / mrchem::orbital::get_electron_number(Phi);
    mrcpp::ComplexFunction pot = calcPotential(rho, poisson_op, abs_prec);
    
    int numAtoms = mol.getNNuclei();
    int numOrbitals = Phi.size();

    Eigen::MatrixXd forces = Eigen::MatrixXd::Zero(numAtoms, 3);
    Vector3d center;
    std::array<double, 3> coord;

    std::filesystem::path p = __FILE__;
    std::filesystem::path parent_dir = p.parent_path();
    std::string filename = parent_dir.string() + "/lebvedev.txt";

    double radius = 0.8;
    VectorXd radii(8);
    radii << .5, .55, 0.6, 65, 0.7, .75, 0.8, .85;

    // loop over radii:
    for (int i = 0; i < radii.size(); i++) {
        radius = radii(i);

        for (int iAtom = 0; iAtom < numAtoms; iAtom++) {
            coord = mol.getNuclei()[iAtom].getCoord();
            center << coord[0], coord[1], coord[2];
            LebedevIntegrator integrator(filename, radius, center);
            MatrixXd gridPos = integrator.getPoints();
            VectorXd weights = integrator.getWeights();
            MatrixXd normals = integrator.getNormals();
            std::vector<Matrix3d> xStres = xcStress(mol, rho, gridPos, json_fock, prec);
            std::vector<Matrix3d> kstress = kineticStress(mol, Phi, prec, gridPos);
            std::vector<Matrix3d> mstress = maxwellStress(mol, pot, nabla, gridPos);
            std::vector<Matrix3d> stress(integrator.n);
            for (int i = 0; i < integrator.n; i++){
                stress[i] = xStres[i] + kstress[i] + mstress[i];
                forces.row(iAtom) -= stress[i] * normals.row(i).transpose() * weights(i);
            }
        }
    }
    forces = forces / radii.size();

    for (int iAtom = 0; iAtom < numAtoms; iAtom++) {
        std::cerr << "forces " << forces(iAtom, 0) << " " << forces(iAtom, 1) << " " << forces(iAtom, 2) << std::endl;
    }
    

    nabla.clear();
    return forces;
}

} // namespace surface_force
