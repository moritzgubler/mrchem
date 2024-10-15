/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include <MRCPP/MWOperators>
#include <MRCPP/Parallel>
#include <MRCPP/Printer>
#include <MRCPP/Timer>
#include <MRCPP/utils/details.h>

#include <string>
#include <vector>
#include "utils/PolyInterpolator.h"

#include "core.h"
#include "gto.h"
#include "sad.h"

#include "chemistry/Nucleus.h"

#include "utils/print_utils.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"

#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/two_electron/CoulombOperator.h"
#include "qmoperators/two_electron/XCOperator.h"
#include "qmoperators/one_electron/AZoraPotential.h"
#include "qmoperators/one_electron/ZoraOperator.h"

#include "mrdft/Factory.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

namespace initial_guess {
namespace sad {

void project_atomic_densities(double prec, Density &rho_tot, const Nuclei &nucs, double screen = -1.0);
void project_atomic_densities_new(double prec, Density &rho_tot, const Nuclei &nucs);

} // namespace sad
} // namespace initial_guess

bool initial_guess::sad::setup(OrbitalVector &Phi, double prec, double screen, const Nuclei &nucs, int zeta) {
    if (Phi.size() == 0) return false;

    auto restricted = (orbital::size_singly(Phi)) ? false : true;
    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation ", "Compute initial orbitals");
    print_utils::text(0, "Method      ", "Diagonalize SAD Hamiltonian");
    print_utils::text(0, "Precision   ", print_utils::dbl_to_str(prec, 5, true));
    print_utils::text(0, "Screening   ", print_utils::dbl_to_str(screen, 5, true) + " StdDev");
    print_utils::text(0, "Restricted  ", (restricted) ? "True" : "False");
    print_utils::text(0, "Functional  ", "LDA (SVWN5)");
    print_utils::text(0, "AO basis    ", "Hydrogenic orbitals");
    print_utils::text(0, "Zeta quality", std::to_string(zeta));
    mrcpp::print::separator(0, '~', 2);

    // Make Fock operator contributions
    Timer t_tot, t_lap;
    auto P_p = std::make_shared<mrcpp::PoissonOperator>(*MRA, prec);
    auto D_p = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.0, 0.0);

    mrdft::Factory xc_factory(*MRA);
    xc_factory.setSpin(false);
    xc_factory.setFunctional("SLATERX", 1.0);
    xc_factory.setFunctional("VWN5C", 1.0);
    auto mrdft_p = xc_factory.build();

    MomentumOperator p(D_p);
    NuclearOperator V_nuc(nucs, prec);
    CoulombOperator J(P_p);
    XCOperator XC(mrdft_p);
    RankZeroOperator V = V_nuc + J + XC;

    auto plevel = Printer::getPrintLevel();
    if (plevel == 1) mrcpp::print::header(1, "SAD Initial Guess");
    if (plevel == 1) mrcpp::print::time(1, "Initializing operators", t_lap);

    // Compute Coulomb density
    t_lap.start();
    Density &rho_j = J.getDensity();
    // initial_guess::sad::project_atomic_densities(prec, rho_j, nucs, screen);
    initial_guess::sad::project_atomic_densities_new(prec, rho_j, nucs);

    // Compute XC density
    Density &rho_xc = XC.getDensity(DensityType::Total);
    mrcpp::cplxfunc::deep_copy(rho_xc, rho_j);
    if (plevel == 1) mrcpp::print::time(1, "Projecting GTO density", t_lap);

    // Project AO basis of hydrogen functions
    t_lap.start();
    OrbitalVector Psi;
    initial_guess::core::project_ao(Psi, prec, nucs, zeta);
    if (plevel == 1) mrcpp::print::time(1, "Projecting Hydrogen AOs", t_lap);

    if (plevel == 2) mrcpp::print::header(2, "Building Fock operator");
    t_lap.start();
    p.setup(prec);
    V.setup(prec);
    if (plevel == 2) mrcpp::print::footer(2, t_lap, 2);
    if (plevel == 1) mrcpp::print::time(1, "Building Fock operator", t_lap);

    // Compute Fock matrix
    mrcpp::print::header(2, "Diagonalizing Fock matrix");
    ComplexMatrix U = initial_guess::core::diagonalize(Psi, p, V);

    // Rotate orbitals and fill electrons by Aufbau
    t_lap.start();
    auto Phi_a = orbital::disjoin(Phi, SPIN::Alpha);
    auto Phi_b = orbital::disjoin(Phi, SPIN::Beta);
    initial_guess::core::rotate_orbitals(Phi, prec, U, Psi);
    initial_guess::core::rotate_orbitals(Phi_a, prec, U, Psi);
    initial_guess::core::rotate_orbitals(Phi_b, prec, U, Psi);
    Phi = orbital::adjoin(Phi, Phi_a);
    Phi = orbital::adjoin(Phi, Phi_b);
    p.clear();
    V.clear();

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::footer(1, t_tot, 2);
    return true;
}

bool initial_guess::sad::setup(OrbitalVector &Phi, double prec, double screen, const Nuclei &nucs) {
    if (Phi.size() == 0) return false;
    bool zora = true;

    auto restricted = (orbital::size_singly(Phi)) ? false : true;
    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation ", "Compute initial orbitals");
    print_utils::text(0, "Method      ", "Diagonalize SAD Hamiltonian");
    print_utils::text(0, "Precision   ", print_utils::dbl_to_str(prec, 5, true));
    print_utils::text(0, "Screening   ", print_utils::dbl_to_str(screen, 5, true) + " StdDev");
    print_utils::text(0, "Restricted  ", (restricted) ? "True" : "False");
    print_utils::text(0, "Functional  ", "LDA (SVWN5)");
    print_utils::text(0, "AO basis    ", "3-21G");
    mrcpp::print::separator(0, '~', 2);

    // Make Fock operator contributions
    Timer t_tot, t_lap;
    auto P_p = std::make_shared<mrcpp::PoissonOperator>(*MRA, prec);
    auto D_p = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.0, 0.0);

    mrdft::Factory xc_factory(*MRA);
    xc_factory.setSpin(false);
    xc_factory.setFunctional("SLATERX", 1.0);
    xc_factory.setFunctional("VWN5C", 1.0);
    auto mrdft_p = xc_factory.build();

    int adap = 0;
    double c = 137.035999084;
    AZoraPotential azoraPot(nucs, adap, prec, getAzoraDir(), false);
    ZoraOperator zoraChi(azoraPot, c, prec, false);
    zoraChi.setup(prec);

    MomentumOperator p(D_p);
    NuclearOperator V_nuc(nucs, prec);
    CoulombOperator J(P_p);
    XCOperator XC(mrdft_p);
    RankZeroOperator V = V_nuc + J + XC;

    auto plevel = Printer::getPrintLevel();
    if (plevel == 1) mrcpp::print::header(1, "SAD Initial Guess");
    if (plevel == 1) mrcpp::print::time(1, "Initializing operators", t_lap);

    // Compute Coulomb density
    t_lap.start();
    Density &rho_j = J.getDensity();
    // initial_guess::sad::project_atomic_densities(prec, rho_j, nucs, screen);
    initial_guess::sad::project_atomic_densities_new(prec, rho_j, nucs);

    // Compute XC density
    Density &rho_xc = XC.getDensity(DensityType::Total);
    mrcpp::cplxfunc::deep_copy(rho_xc, rho_j);
    if (plevel == 1) mrcpp::print::time(1, "Projecting GTO density", t_lap);

    // Project AO basis of hydrogen functions
    t_lap.start();
    OrbitalVector Psi;
    initial_guess::gto::project_ao(Psi, prec, nucs);
    
    if (plevel == 1) mrcpp::print::time(1, "Projecting GTO AOs", t_lap);
    if (plevel == 2) mrcpp::print::header(2, "Building Fock operator");
    t_lap.start();
    p.setup(prec);
    V.setup(prec);
    if (plevel == 2) mrcpp::print::footer(2, t_lap, 2);
    if (plevel == 1) mrcpp::print::time(1, "Building Fock operator", t_lap);

    // Compute Fock matrix
    mrcpp::print::header(2, "Diagonalizing Fock matrix");

    ComplexMatrix U;
    if (zora) {
        U = initial_guess::core::diagonalize_with_zora(Psi, p, V, zoraChi);
    } else {
        U = initial_guess::core::diagonalize(Psi, p, V);
    }

    // Rotate orbitals and fill electrons by Aufbau
    t_lap.start();
    auto Phi_a = orbital::disjoin(Phi, SPIN::Alpha);
    auto Phi_b = orbital::disjoin(Phi, SPIN::Beta);
    initial_guess::core::rotate_orbitals(Phi, prec, U, Psi);
    initial_guess::core::rotate_orbitals(Phi_a, prec, U, Psi);
    initial_guess::core::rotate_orbitals(Phi_b, prec, U, Psi);
    Phi = orbital::adjoin(Phi, Phi_a);
    Phi = orbital::adjoin(Phi, Phi_b);
    p.clear();
    V.clear();
    zoraChi.clear();

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::footer(1, t_tot, 2);
    return true;
}

void initial_guess::sad::project_atomic_densities_new(double prec, Density &rho, const Nuclei &nucs){
    std::string source_dir = HIRSHFELD_SOURCE_DIR;
    std::string install_dir = HIRSHFELD_INSTALL_DIR;
    std::string data_dir = "";
    // check if data_dir exists
    if (std::filesystem::exists(install_dir)) {
        data_dir = install_dir + "/lda/";
    } else if (std::filesystem::exists(source_dir)) {
        data_dir = source_dir + "/lda/";
    } else {
        MSG_ABORT("Hirshfeld data directory not found");
    }
    std::vector<interpolation_utils::PolyInterpolator> atomic_densities;

    for (const auto &nuc : nucs) {
        std::string element = nuc.getElement().getSymbol();
        std::string file = data_dir + "/" + element + ".density";

        Eigen::VectorXd rhoGrid, rGrid;
        mrchem::density::readAtomicDensity(file, rGrid, rhoGrid);

        interpolation_utils::PolyInterpolator atomic_density(rGrid, rhoGrid);
        atomic_densities.push_back(atomic_density);
    }
    auto rho_analytic = [atomic_densities, nucs](const mrcpp::Coord<3> &r) {
        double rho = 0.0;
        for (int i = 0; i < nucs.size(); i++) {
            mrcpp::Coord<3> nucPos = nucs[i].getCoord();
            double rr = std::sqrt((r[0] - nucPos[0]) * (r[0] - nucPos[0])
                + (r[1] - nucPos[1]) * (r[1] - nucPos[1])
                + (r[2] - nucPos[2]) * (r[2] - nucPos[2]));
            rho += atomic_densities[i].evalfLeftNoRightZero(rr);
        }
        return rho;
    };
    
    mrcpp::ComplexFunction rho_MW;
    mrcpp::cplxfunc::project(rho_MW, rho_analytic, mrcpp::NUMBER::Real, prec);

    mrcpp::ComplexDouble integral = rho_MW.integrate();
    double norm = integral.real();
    std::cout << "Integral of initial charge density: " << norm << std::endl;
    rho.add(1.0, rho_MW);
    
}

void initial_guess::sad::project_atomic_densities(double prec, Density &rho_tot, const Nuclei &nucs, double screen) {
    auto pprec = Printer::getPrecision();
    auto w0 = Printer::getWidth() - 1;
    auto w1 = 5;
    auto w2 = 8;
    auto w3 = w0 / 3;
    auto w4 = w0 - (w1 + w2 + 2 * w3);

    std::stringstream o_head;
    o_head << std::setw(w1) << "N";
    o_head << std::setw(w2) << "Atom";
    o_head << std::setw(w4) << " ";
    o_head << std::setw(w3) << "Nuclear charge";
    o_head << std::setw(w3) << "Electron charge";

    mrcpp::print::header(2, "Projecting GTO density");
    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    auto crop_prec = (mrcpp::mpi::numerically_exact) ? -1.0 : prec;
    std::string sad_path;
    for (auto n : {sad_basis_source_dir(), sad_basis_install_dir()}) {
        auto trimmed = print_utils::rtrim_copy(n);
        if (mrcpp::details::directory_exists(trimmed)) {
            sad_path = trimmed;
            break;
        }
    }

    Timer t_tot;
    Density rho_loc(false);
    rho_loc.alloc(NUMBER::Real);
    rho_loc.real().setZero();

    Timer t_loc;
    auto N_nucs = nucs.size();
    DoubleVector charges = DoubleVector::Zero(2 * N_nucs);
    for (int k = 0; k < N_nucs; k++) {
        if (mrcpp::mpi::wrk_rank != k % mrcpp::mpi::wrk_size) continue;

        const std::string &sym = nucs[k].getElement().getSymbol();
        std::stringstream o_bas, o_dens;
        o_bas << sad_path << "/" << sym << ".bas";
        o_dens << sad_path << "/" << sym << ".dens";

        Density rho_k = initial_guess::gto::project_density(prec, nucs[k], o_bas.str(), o_dens.str(), screen);
        rho_loc.add(1.0, rho_k);
        rho_loc.crop(crop_prec);

        charges[k] = nucs[k].getCharge();
        charges[N_nucs + k] = rho_k.integrate().real();
    }
    t_loc.stop();
    Timer t_com;
    mrcpp::mpi::allreduce_vector(charges, mrcpp::mpi::comm_wrk);
    density::allreduce_density(prec, rho_tot, rho_loc);
    t_com.stop();

    for (int k = 0; k < N_nucs; k++) {
        std::stringstream o_row;
        o_row << std::setw(w1) << k;
        o_row << std::setw(w2) << nucs[k].getElement().getSymbol();
        o_row << std::setw(w4) << " ";
        o_row << std::setw(w3) << print_utils::dbl_to_str(charges[k], 2 * pprec, false);
        o_row << std::setw(w3) << print_utils::dbl_to_str(charges[N_nucs + k], 2 * pprec, false);
        println(2, o_row.str());
    }

    auto tot_nuc = charges.segment(0, N_nucs).sum();
    auto tot_rho = charges.segment(N_nucs, N_nucs).sum();

    std::stringstream o_row;
    o_row << " Total charge";
    o_row << std::string(w1 + w2 + w4 - 13, ' ');
    o_row << std::setw(w3) << print_utils::dbl_to_str(tot_nuc, 2 * pprec, false);
    o_row << std::setw(w3) << print_utils::dbl_to_str(tot_rho, 2 * pprec, false);

    mrcpp::print::separator(2, '-');
    println(2, o_row.str());
    mrcpp::print::separator(2, '-');
    print_utils::qmfunction(2, "Local density", rho_loc, t_loc);
    print_utils::qmfunction(2, "Allreduce density", rho_tot, t_com);
    mrcpp::print::footer(2, t_tot, 2);
}

} // namespace mrchem
