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
     * @param adap Adaptive parameter from base class.
     * @param prec Precision parameter from base class.
     * @param shared Determines if the base potential is shared.
     */
    AZoraPotential(const Molecule& molecule, int adap, double prec, bool inverse, bool shared = false) 
        : QMPotential(adap, shared), mol(molecule), prec(prec), inverse(inverse) {
        initAzoraPotential(other.inverse);
    }

    /**
     * Copy constructor.
     * @param other The other instance to copy from.
     */
    AZoraPotential(const AZoraPotential& other) 
        : QMPotential(other), mol(other.mol), prec(other.prec), inverse(other.inverse) {
        initAzoraPotential(other.inverse);
    }

    /**
     * Destructor.
     */
    virtual ~AZoraPotential() override = default;

    // Delete copy assignment to prevent copying
    AZoraPotential& operator=(const AZoraPotential&) = delete;

protected:
    Molecule mol; // The molecule representing the potential
    double prec; // The precision parameter
    bool inverse; // Determines if the potential is inverted or not

    /**
     * Initialize the azora potential based on the molecule.
     * This method would typically setup the real and imaginary function trees
     * representing the potential.
     * @param inverse Determines if the kapppa function is inverted or not.
     */
    void initAzoraPotential(bool inverse = false) {

        int n = mol.getNNuclei();
        Nuclei nucs = mol.getNuclei();
        Eigen::VectorXd rGrid, vZora, kappa;

        std::vector<RadInterpolater> kappaSplines;


        for (int i = 0; i < n; i++) {
            RadInterpolater kappaSpline(nucs[i].getElement().getSymbol());
            kappaSplines.push_back(kappaSpline);
        }

        // Create lambda function that sums up all the atomic dampening functions
        auto k = [kappaSplines, inverse](const mrcpp::Coord<3>& r) {
            double V = 1.0;
            // Loop over all atoms:
            for (int i = 0; i < kappaSplines.size(); i++) {
                double r_i = (r - nucs[i].getPos()).norm();
                V += kappaSplines[i](r_i) - 1.0;
            }
            if (inverse) V = 1.0 / V;
            return V;
        };

        mrcpp::ComplexFunction vtot;
        mrcpp::cplxfunc::project(vtot, V, mrcpp::NUMBER::Real, prec);
        this->add(1.0, vtot);

    }


    // Override the base class methods if needed to adapt them for the AZora potential
};

} // namespace mrchem
