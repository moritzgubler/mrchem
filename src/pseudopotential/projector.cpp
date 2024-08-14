
#include "pseudopotential/projector.h"
#include <math.h>
// #include "mrchem.h"

/**
 * @brief Constructs a Projector object.
 *
 * This constructor initializes a Projector object with the given parameters.
 *
 * @param pos The position vector of the corresponding atom.
 * @param rl The radial length of the projector.
 * @param i Vector index of the projector.
 * @param l The angular momentum quantum number of the projector.
 * @param m The magnetic quantum number of the projector.
 * @param prec The precision of the projector.
 */
ProjectorFunction::ProjectorFunction(mrcpp::Coord<3> pos, double rl, int i, int l, int m, double prec) {
    std::cout << "Constructing ProjectorFunction" << std::endl;
    this->pos = pos;
    this->rl = rl;
    this->i = i;
    this->l = l;
    this->m = m;
    this->prec = prec;
    // select the spherical harmonic function based on the angular momentum and magnetic quantum number
    switch_sperics(l, m);
    double prefactor = std::sqrt(2.0) / (std::pow(rl, l + (4 * i - 2) / 2) * std::sqrt(tgamma( l + (4 * i - 2) / 2 )) );

    auto project_analytic = [this, prefactor](const std::array<double, 3> &r) -> double {
        std::array<double, 3> rprime = {r[0] - this->pos[0], r[1] - this->pos[1], r[2] - this->pos[2]};
        double normr = std::sqrt( rprime[0] * rprime[0] + rprime[1] * rprime[1] + rprime[2] * rprime[2]);
        return prefactor * std::pow(normr, this->l + 2 * (this->i - 1)) * std::exp(- normr * normr / (this->rl * this-> rl) ) * this->s(rprime, normr);
    };
    auto op = (*this);
    mrcpp::cplxfunc::project(op, project_analytic, mrcpp::NUMBER::Real, prec);
    std::cout << "ProjectorFunction constructed in constructor" << std::endl;
}

void ProjectorFunction::switch_sperics(int l, int m){
    switch (l) {
        case 0:
            switch (m) {
                case 0:
                    s = s0;
                    break;
            }
            break;
        case 1:
            switch (m) {
                case -1:
                    s = s1m1;
                    break;
                case 0:
                    s = s10;
                    break;
                case 1:
                    s = s11;
                    break;
            }
            break;
        case 2:
            switch (m) {
                case -2:
                    s = s2m2;
                    break;
                case -1:
                    s = s2m1;
                    break;
                case 0:
                    s = s20;
                    break;
                case 1:
                    s = s21;
                    break;
                case 2:
                    s = s22;
                    break;
            }
            break;
        case 3:
            switch (m) {
                case -3:
                    s = s3m3;
                    break;
                case -2:
                    s = s3m2;
                    break;
                case -1:
                    s = s3m1;
                    break;
                case 0:
                    s = s30;
                    break;
                case 1:
                    s = s31;
                    break;
                case 2:
                    s = s32;
                    break;
                case 3:
                    s = s33;
                    break;
            }
            break;
    }
}