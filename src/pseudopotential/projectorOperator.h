#pragma once

#include "tensor/RankZeroOperator.h"
#include "mrchem.h"
#include "pseudopotential/pseudopotential.h"
#include "pseudopotential/projector.h"
#include "qmfunctions/Orbital.h"
#include "chemistry/Molecule.h"

class Projector: public mrchem::RankZeroOperator {

    std::vector<PseudopotentialData> pp;
    std::vector<Projector> proj;

public:
    Projector(mrchem::Molecule molecule, std::vector<PseudopotentialData> pp){
    }

mrchem::Orbital apply(mrchem::Orbital phi_p) {
    return phi_p;
}

mrchem::Orbital dagger(mrchem::Orbital phi_p) {
    return apply(phi_p);
}

mrchem::QMOperatorVector apply(std::shared_ptr<mrchem::QMOperator> &O) {
    NOT_IMPLEMENTED_ABORT;
}

};