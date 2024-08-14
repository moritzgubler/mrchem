#pragma once

#include "tensor/RankZeroOperator.h"
#include "mrchem.h"
#include "pseudopotential/pseudopotential.h"
#include "pseudopotential/projector.h"
#include "qmfunctions/Orbital.h"
#include "chemistry/Molecule.h"
#include "qmfunctions/qmfunction_utils.h"
#include <string>

class magneticQuantumNumberProjector {
    public:
    std::vector<ProjectorFunction> iProj;
    int nProj;
};

class angularMomentumProjector {
    public:
    std::vector<magneticQuantumNumberProjector> mProj;
    int nM;
};

class AtomProjector {
    public:
    std::vector<angularMomentumProjector> lProj;
    int numberOfAngMom;
};

class ProjectorOperator final : public mrchem::RankZeroOperator {

    std::vector<PseudopotentialData> pp;
    std::vector<AtomProjector> proj;
    double prec;

public:
    ProjectorOperator(mrchem::Molecule &molecule, double prec){

        mrchem::Nuclei nucs = molecule.getNuclei();
        for (int i = 0; i < nucs.size(); i++){
            std::string elem = nucs[i].getElement().getSymbol();
            std::string fname = "psppar." + elem;
            this->pp.push_back(PseudopotentialData(fname));
            pp[i].print();
        }

        std::cout << "Pseudopotential data loaded" << std::endl;

        this->pp = pp;
        this->prec = prec;
        int npp = 0;

        // loop over all atoms and create projectors
        for (int i = 0; i < molecule.getNNuclei(); i++) {
            std::cout << "Creating projectors for atom " << i << std::endl;
            mrcpp::Coord<3> pos = molecule.getNuclei()[i].getCoord();
            proj.push_back(AtomProjector());
            for (int l = 0; l < pp[i].nsep; l++) {
                proj[i].lProj.push_back(angularMomentumProjector());
                for (int m = -l; m <= l; m++) {
                    proj[i].lProj[l].mProj.push_back(magneticQuantumNumberProjector());
                    for (int idim = 0; idim < pp[i].dim_h[l]; idim++){
                        // proj.push_back(ProjectorFunction(pos, pp[i].rl[l], isep, l, m, prec));
                        std::cout << "Creating ProjectorFunction " << l << " " << m << " " << idim << std::endl;
                        proj[i].lProj[l].mProj[m].iProj.push_back(ProjectorFunction(pos, pp[i].rl[l], idim, l, m, prec));
                        std::cout << "ProjectorFunction constructed " << i << std::endl;
                        std::cout << "i = " << i << std::endl;
                        std::cout << "nsep = " << pp[i].nsep << std::endl;
                        proj[i].numberOfAngMom = pp[i].nsep;
                        std::cout << "ProjectorFunction added to projector" << std::endl;
                        proj[i].lProj[l].nM = 2*pp[i].nsep + 1;
                        std::cout << "ProjectorFunction added to projector" << std::endl;
                        proj[i].lProj[l].mProj[m].nProj = pp[i].dim_h[l];
                        std::cout << "End of loop" << std::endl;
                        npp++;
                    }
                }
            }
        }
        std::cout << "ProjectorOperator constructed" << std::endl;
    }

    void setup(double prec) {
        this->prec = prec;
    }

    void clear() {
    }


mrchem::Orbital apply(mrchem::Orbital phi) {
    std::cout << "Applying projector operator" << std::endl;
    // loop over all atoms
    ComplexDouble dotComplex;

    std::vector<ComplexDouble> complexCoefficients;
    mrchem::ComplexFunctionVector complexFunctionVector;

    for (int iat = 0; iat < proj.size(); iat++) {
        // loop over all angular momenta
        for (int l = 0; l < pp[iat].nsep; l++){
            // loop over all magnetic quantum numbers
            for (int m = -l; m <= l; m++){
                // loop over all projectors
                Eigen::VectorXd dot_products(pp[iat].dim_h[l]);
                for (int ip = 0; ip < pp[iat].dim_h[l]; ip++){
                    dotComplex = mrchem::qmfunction::dot(phi, proj[iat].lProj[l].mProj[m].iProj[ip]);
                    dot_products(ip) = dotComplex.real();
                }
                dot_products = pp[iat].h[l] * dot_products;
                // loop over all projectors
                for (int ip = 0; ip < pp[iat].dim_h[l]; ip++){
                    complexCoefficients.push_back(dot_products(ip));
                    complexFunctionVector.push_back(proj[iat].lProj[l].mProj[m].iProj[ip]);
                }
            }
        }
        
    }
    // convert complexCoefficients to Eigen Vector:
    mrchem::ComplexVector complexCoefficientsEigen = Eigen::Map<Eigen::VectorXcd>(complexCoefficients.data(), complexCoefficients.size());

    mrchem::Orbital result;
    mrchem::qmfunction::linear_combination(result, complexCoefficientsEigen, complexFunctionVector, prec);
    return result;
}

mrchem::Orbital dagger(mrchem::Orbital phi) {
    return apply(phi);
}

mrchem::QMOperatorVector apply(std::shared_ptr<mrchem::QMOperator> &O) {
    NOT_IMPLEMENTED_ABORT;
}

};