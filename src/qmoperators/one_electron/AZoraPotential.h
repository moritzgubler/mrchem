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
    AZoraPotential(Nuclei nucs, int adap, double prec, bool inverse, std::string mode, bool shared = false) 
        : QMPotential(adap, shared) {
        this->nucs = nucs;
        this->prec = prec;
        this->inverse = inverse;
        this->mode = mode;
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
        this->inverse = other.inverse;
        this->mode = other.mode;
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
    bool inverse; // Whether to use the inverse of the damping function
    std::string mode; // The mode of the potential, either 'kappa' or 'potential'

    /**
     * Initialize the azora potential based on the molecule.
     * This method would typically setup the real and imaginary function trees
     * representing the potential.
     */
    void initAzoraPotential() {

        int n = nucs.size();
        Eigen::VectorXd rGrid, vZora, kappa;

        std::vector<RadInterpolater> dampRel;


        for (int i = 0; i < n; i++) {
            RadInterpolater kappaSpline(nucs[i].getElement().getSymbol(), this->mode);
            dampRel.push_back(kappaSpline);
        }

        std::cerr << " Sizs of dampRel: " << dampRel.size() << "\n";

        mrcpp::ComplexFunction vtot;
        if (mode == "kappa") {
        auto k = [dampRel, this](const mrcpp::Coord<3>& r) {
            double V = 0.0;
            // Loop over all atoms:
            for (int i = 0; i < dampRel.size(); i++) {
                mrcpp::Coord<3> r_i = nucs[i].getCoord();
                double rr = std::sqrt((r[0] - r_i[0]) * (r[0] - r_i[0]) + (r[1] - r_i[1]) * (r[1] - r_i[1]) + (r[2] - r_i[2]) * (r[2] - r_i[2]));
                V += dampRel[i].evalf(rr) - .5;
            }
            V += .5;
            if (this->inverse) V = 1.0 / V;
            return V;
        };
        mrcpp::cplxfunc::project(vtot, k, mrcpp::NUMBER::Real, prec);
        } else if (mode == "potential") {
        auto k = [dampRel, this](const mrcpp::Coord<3>& r) {
            double V = 0.0;
            // Loop over all atoms:
            for (int i = 0; i < dampRel.size(); i++) {
                mrcpp::Coord<3> r_i = nucs[i].getCoord();
                double rr = std::sqrt((r[0] - r_i[0]) * (r[0] - r_i[0]) + (r[1] - r_i[1]) * (r[1] - r_i[1]) + (r[2] - r_i[2]) * (r[2] - r_i[2]));
                V += dampRel[i].evalf(rr);
            }
            return V;
        };
        mrcpp::cplxfunc::project(vtot, k, mrcpp::NUMBER::Real, prec);
        }
        this->add(1.0, vtot);

    }
};

} // namespace mrchem
