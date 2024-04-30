#pragma once

#include "tensor/RankZeroOperator.h"
#include "mrchem.h"
#include "pseudopotential/pseudopotential.h"
#include "qmfunctions/Orbital.h"
#include "chemistry/Molecule.h"

class Projector: public mrchem::RankZeroOperator {

public:
    Projector(mrchem::Molecule molecule, PseudopotentialData pp){
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