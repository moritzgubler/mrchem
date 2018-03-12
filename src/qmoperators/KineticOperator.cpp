#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "KineticOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

ComplexDouble KineticOperator::operator()(Orbital bra, Orbital ket) {
    RankZeroTensorOperator &p_x = this->p[0];
    RankZeroTensorOperator &p_y = this->p[1];
    RankZeroTensorOperator &p_z = this->p[2];

    Orbital bra_x = p_x(bra);
    Orbital ket_x = p_x(ket);
    ComplexDouble T_x = orbital::dot(bra_x, ket_x);
    bra_x.free();
    ket_x.free();

    Orbital bra_y = p_y(bra);
    Orbital ket_y = p_y(ket);
    ComplexDouble T_y = orbital::dot(bra_y, ket_y);
    bra_y.free();
    ket_y.free();

    Orbital bra_z = p_z(bra);
    Orbital ket_z = p_z(ket);
    ComplexDouble T_z = orbital::dot(bra_z, ket_z);
    bra_z.free();
    ket_z.free();

    return 0.5*(T_x + T_y + T_z);
}

ComplexDouble KineticOperator::dagger(Orbital bra, Orbital ket) {
    NOT_IMPLEMENTED_ABORT;
}

ComplexMatrix KineticOperator::operator()(OrbitalVector &bra, OrbitalVector &ket) {
    Timer timer;
    Printer::printHeader(1, "Compute Kinetic Matrix Elements");

    RankZeroTensorOperator &p_x = this->p[0];
    RankZeroTensorOperator &p_y = this->p[1];
    RankZeroTensorOperator &p_z = this->p[2];

    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix T_x = ComplexMatrix::Zero(Ni, Nj);
    ComplexMatrix T_y = ComplexMatrix::Zero(Ni, Nj);
    ComplexMatrix T_z = ComplexMatrix::Zero(Ni, Nj);
    {
        Timer timer;
        OrbitalVector dKet = p_x(ket);
        for (int i = 0; i < Ni; i++) {
            Orbital dBra_i = p_x(bra[i]);
            for (int j = 0; j < Nj; j++) {
                T_x(i,j) = orbital::dot(dBra_i, dKet[j]);
            }
            dBra_i.free();
        }
        orbital::free(dKet);
        timer.stop();
        double t = timer.getWallTime();
        Printer::printDouble(1, "T_x", t);
    }
    {
        Timer timer;
        OrbitalVector dKet = p_y(ket);
        for (int i = 0; i < Ni; i++) {
            Orbital dBra_i = p_y(bra[i]);
            for (int j = 0; j < Nj; j++) {
                T_x(i,j) = orbital::dot(dBra_i, dKet[j]);
            }
            dBra_i.free();
        }
        orbital::free(dKet);
        timer.stop();
        double t = timer.getWallTime();
        Printer::printDouble(1, "T_y", t);
    }
    {
        Timer timer;
        OrbitalVector dKet = p_z(ket);
        for (int i = 0; i < Ni; i++) {
            Orbital dBra_i = p_z(bra[i]);
            for (int j = 0; j < Nj; j++) {
                T_x(i,j) = orbital::dot(dBra_i, dKet[j]);
            }
            dBra_i.free();
        }
        orbital::free(dKet);
        timer.stop();
        double t = timer.getWallTime();
        Printer::printDouble(1, "T_z", t);
    }
    timer.stop();
    Printer::printFooter(1, timer, 2);

    return 0.5*(T_x + T_y + T_z);
}

ComplexMatrix KineticOperator::dagger(OrbitalVector &bra, OrbitalVector &ket) {
    NOT_IMPLEMENTED_ABORT;
}

} //namespace mrchem