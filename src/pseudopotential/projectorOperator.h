#pragma once

#include "tensor/RankZeroOperator.h"
#include "mrchem.h"
#include "pseudopotential/pseudopotential.h"
#include "pseudopotential/projector.h"
#include "qmfunctions/Orbital.h"
#include "chemistry/Molecule.h"
#include "qmfunctions/qmfunction_utils.h"

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

class ProjectorOperator: public mrchem::RankZeroOperator {

    std::vector<PseudopotentialData> pp;
    std::vector<AtomProjector> proj;
    double prec;

public:
    ProjectorOperator(mrchem::Molecule molecule, std::vector<PseudopotentialData> pp, double prec){

        this->pp = pp;
        this->prec = prec;
        int npp = 0;

        // loop over all atoms and create projectors
        for (int i = 0; i < molecule.getNNuclei(); i++) {
            mrcpp::Coord<3> pos = molecule.getNuclei()[i].getCoord();
            proj.push_back(AtomProjector());
            for (int l = 0; l < pp[i].nsep; l++) {
                proj[i].lProj.push_back(angularMomentumProjector());
                for (int m = -l; m <= l; m++) {
                    proj[i].lProj[l].mProj.push_back(magneticQuantumNumberProjector());
                    for (int idim = 0; idim < pp[i].dim_h[l]; idim++){
                        // proj.push_back(ProjectorFunction(pos, pp[i].rl[l], isep, l, m, prec));
                        proj[i].lProj[l].mProj[m].iProj.push_back(ProjectorFunction(pos, pp[i].rl[l], idim, l, m, prec));
                        proj[i].numberOfAngMom = pp[i].nsep;
                        proj[i].lProj[l].nM = 2*pp[i].nsep + 1;
                        proj[i].lProj[l].mProj[m].nProj = pp[i].dim_h[l];
                        npp++;
                    }
                }
            }
        }

    }

mrchem::Orbital apply(mrchem::Orbital phi) {
    // loop over all atoms
    for (int iat = 0; iat < proj.size(); iat++) {
        // loop over all angular momenta
    }
    
    return phi;
}

mrchem::Orbital dagger(mrchem::Orbital phi) {
    return apply(phi);
}

mrchem::QMOperatorVector apply(std::shared_ptr<mrchem::QMOperator> &O) {
    NOT_IMPLEMENTED_ABORT;
}

};