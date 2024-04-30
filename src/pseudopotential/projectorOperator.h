#pragma once

#include "tensor/RankZeroOperator.h"
#include "mrchem.h"
#include "pseudopotential/pseudopotential.h"
#include "pseudopotential/projector.h"
#include "qmfunctions/Orbital.h"
#include "chemistry/Molecule.h"

class ProjectorOperator: public mrchem::RankZeroOperator {

    std::vector<PseudopotentialData> pp;
    std::vector<ProjectorFunction> proj;
    double prec;

public:
    ProjectorOperator(mrchem::Molecule molecule, std::vector<PseudopotentialData> pp, double prec){

        this->pp = pp;
        this->prec = prec;

        // loop over all atoms and create projectors
        for (int i = 0; i < molecule.getNNuclei(); i++) {
            mrcpp::Coord<3> pos = molecule.getNuclei()[i].getCoord();
            for (int l = 0; l < pp[i].nProjectors; l++) {
                for (int m = -l; m <= l; m++) {
                    for (int isep = 0; isep < pp[i].nsep; isep++){
                        proj.push_back(ProjectorFunction(pos, pp[i].rl[l], isep, l, m, prec));
                    }
                    
                }
            }
        }

    }

// mrchem::Orbital apply(mrchem::Orbital phi_p) {
//     return phi_p;
// }

// mrchem::Orbital dagger(mrchem::Orbital phi_p) {
//     return apply(phi_p);
// }

// mrchem::QMOperatorVector apply(std::shared_ptr<mrchem::QMOperator> &O) {
//     NOT_IMPLEMENTED_ABORT;
// }

};