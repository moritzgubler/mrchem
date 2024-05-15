#pragma once

#include "mrchem.h"
#include "chemistry/Molecule.h"
#include "chemistry/Nucleus.h"
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <string>
#include "qmoperators/one_electron/AZora/radialInterpolater.h"


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
    AZoraPotential(const std::shared_ptr<Molecule> molecule, int adap, double prec, bool shared = false) 
        : QMPotential(adap, shared) {
        this->mol = molecule;
        this->prec = prec;
        initAzoraPotential();
    }

    /**
     * Copy constructor.
     * @param other The other instance to copy from.
     */
    AZoraPotential(const AZoraPotential& other) 
        : QMPotential(other) {
        this->mol = other.mol;
        this->prec = other.prec;
        initAzoraPotential();
    }

    /**
     * Destructor.
     */
    virtual ~AZoraPotential() override = default;

    // Delete copy assignment to prevent copying
    AZoraPotential& operator=(const AZoraPotential&) = delete;

protected:
    std::shared_ptr<const Molecule> mol; // The molecule representing the potential
    double prec; // The precision parameter

    /**
     * Initialize the azora potential based on the molecule.
     * This method would typically setup the real and imaginary function trees
     * representing the potential.
     */
    void initAzoraPotential() {

        int n = mol->getNNuclei();
        Nuclei nucs = mol->getNuclei();
        Eigen::VectorXd rGrid, vZora, kappa;

        std::vector<RadInterpolater> aZoraPotentialSplines;


        for (int i = 0; i < n; i++) {
            RadInterpolater kappaSpline(nucs[i].getElement().getSymbol());
            aZoraPotentialSplines.push_back(kappaSpline);
        }

        // Create lambda function that sums up all the atomic dampening functions
        auto k = [aZoraPotentialSplines, nucs](const mrcpp::Coord<3>& r) {
            double V = 1.0;
            // Loop over all atoms:
            for (int i = 0; i < aZoraPotentialSplines.size(); i++) {
                mrcpp::Coord<3> r_i = nucs[i].getCoord();
                double rr = std::sqrt((r[0] - r_i[0]) * (r[0] - r_i[0]) + (r[1] - r_i[1]) * (r[1] - r_i[1]) + (r[2] - r_i[2]) * (r[2] - r_i[2]));
                V += aZoraPotentialSplines[i].evalf(rr) - 1.0;
            }
            return V;
        };

        mrcpp::ComplexFunction vtot;
        mrcpp::cplxfunc::project(vtot, k, mrcpp::NUMBER::Real, prec);
        this->add(1.0, vtot);

    }


    // Override the base class methods if needed to adapt them for the AZora potential
};

} // namespace mrchem
