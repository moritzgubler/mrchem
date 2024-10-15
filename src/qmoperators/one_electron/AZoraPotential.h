#pragma once

#include "mrchem.h"
#include "chemistry/Molecule.h"
#include "chemistry/Nucleus.h"
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <string>
#include "qmoperators/one_electron/AZora/RadialInterpolater.h"


namespace mrchem {

/** @class AZora
 *
 * @brief Operator defining the AZora potential based on molecular data.
 *
 * Inherits from QMPotential and adds functionality to utilize an mrchem::Molecule
 * for constructing the potential.
 *
 */
class AZoraPotential : public QMPotential {
public:
    /**
     * Constructor that takes a molecule and initializes the azora potential.
     * @param molecule The molecule used to construct the potential.
     * @param adap Adaptive parameter from QMPotential.
     * @param prec Precision parameter from base class.
     * @param shared Determines if the base potential is shared.
     */
    AZoraPotential(Nuclei nucs, int adap, double prec, std::string azora_dir, bool shared = false, double c = 137.035999084) 
        : QMPotential(adap, shared) {
        this->nucs = nucs;
        this->prec = prec;
        this->c = c;
        this->azora_dir = azora_dir;
        initAzoraPotential();
    }

    /**
     * Copy constructor.
     * @param other The other instance to copy from.
     */
    AZoraPotential(const AZoraPotential& other) 
        : QMPotential(other) {
        this->nucs = other.nucs;
        this->prec = other.prec;
        this->c = other.c;
        this->azora_dir = other.azora_dir;
        initAzoraPotential();
    }

    /**
     * Destructor.
     */
    virtual ~AZoraPotential() override = default;

    // Delete copy assignment to prevent copying
    AZoraPotential& operator=(const AZoraPotential&) = delete;

protected:
    Nuclei nucs; // The nuclei of the molecule
    double prec; // The precision parameter
    double c;    // The speed of light
    std::string azora_dir; // The directory containing the azora potential data

    /**
     * Initialize the azora potential based on the molecule.
     * This method would typically setup the real and imaginary function trees
     * representing the potential.
     */
    void initAzoraPotential() {

        int n = nucs.size();
        Eigen::VectorXd rGrid, vZora, kappa;

        std::vector<RadInterpolater> atomicPotentials;

        std::string mode = "potential";
        for (int i = 0; i < n; i++) {
            RadInterpolater potentialSpline(nucs[i].getElement().getSymbol(), this->azora_dir, mode);
            atomicPotentials.push_back(potentialSpline);
        }


        mrcpp::ComplexFunction vtot;
        auto kappa_analytic = [atomicPotentials, this](const mrcpp::Coord<3>& r) {
            double V = 0.0;
            // Loop over all atoms:
            for (int i = 0; i < atomicPotentials.size(); i++) {
                mrcpp::Coord<3> r_i = nucs[i].getCoord();
                double rr = std::sqrt((r[0] - r_i[0]) * (r[0] - r_i[0]) + (r[1] - r_i[1]) * (r[1] - r_i[1]) + (r[2] - r_i[2]) * (r[2] - r_i[2]));
                V += atomicPotentials[i].evalf(rr);
            }
            return 1 / (1 - V / (2.0 * c * c)) - 1;
        };
        mrcpp::cplxfunc::project(vtot, kappa_analytic, mrcpp::NUMBER::Real, prec);
        // auto ttt = [](const mrcpp::Coord<3>& r) {
        //     return 1.0;
        // };
        // mrcpp::ComplexFunction const_func;
        // mrcpp::cplxfunc::project(const_func, ttt, mrcpp::NUMBER::Real, prec);
        mrcpp::ComplexDouble one = 1.0;
        // vtot.add(one, const_func);
        this->add(one, vtot);
    }
};

/**
 * Get the directory containing the azora potential data.
 */
inline std::string getAzoraDir() {
    std::string azora_dir_src = AZORA_POTENTIALS_SOURCE_DIR;
    std::string azora_dir_install = AZORA_POTENTIALS_INSTALL_DIR;

    std::string azora_dir_final;
    // check if azora_dir_install directory exists using cpp standard library
    if (std::filesystem::exists(azora_dir_install)) {
        azora_dir_final = azora_dir_install;
    } else {
        // check if azora_dir_src directory exists using cpp standard library
        if (std::filesystem::exists(azora_dir_src)) {
            azora_dir_final = azora_dir_src;
        } else {
            MSG_ABORT("AZORA: No directory provided and no default directories found.");
        }
    }
    return azora_dir_final;
}

} // namespace mrchem