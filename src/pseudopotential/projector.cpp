
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
    int ii = i + 1;
    this->l = l;
    this->m = m;
    this->prec = prec;
    // select the spherical harmonic function based on the angular momentum and magnetic quantum number
    switch_sperics(l, m);
    double prefactor = std::sqrt(2.0) / (std::pow(rl, l + (4.0 * ii - 1) / 2.0) * std::sqrt(tgamma( l + (4.0 * ii - 1.0) / 2.0 )) );
    std::cout << "prefactor: " << prefactor << std::endl;
    std::cout << "rl: " << rl << std::endl;
    std::cout << "pow " << std::pow(rl, l + (4.0 * ii - 1) / 2.0) << std::endl;

    // debug with prints why prefactor is nan:
    std::cout << "prefactor: " << prefactor << std::endl;
    std::cout << "rl: " << rl << std::endl;
    std::cout << "l: " << l << std::endl;
    std::cout << "i: " << ii << std::endl;
    std::cout << "tgamma: " << tgamma( l + (4.0 * ii - 1.0) / 2.0 ) << std::endl;
    std::cout << "arg of tgamma: " << l + (4.0 * ii - 1.0) / 2.0 << std::endl;



    auto project_analytic = [this, prefactor, ii](const std::array<double, 3> &r) -> double {
        std::array<double, 3> rprime = {r[0] - this->pos[0], r[1] - this->pos[1], r[2] - this->pos[2]};
        double normr = std::sqrt( rprime[0] * rprime[0] + rprime[1] * rprime[1] + rprime[2] * rprime[2]);
        return prefactor * std::pow(normr, this->l + 2 * (ii - 1)) * std::exp(- 0.5 * normr * normr / (this->rl * this-> rl) ) * this->s(rprime, normr);
    };
    auto op = (*this);
    mrcpp::cplxfunc::project(op, project_analytic, mrcpp::NUMBER::Real, prec);

    mrcpp::Coord<3> r = {0.1, 0.1, 0.1};
    std::cout << "ProjectorFunction at origin: " << this->real().evalf(r) << std::endl;
    std::cout << "analytic at origin: " << project_analytic(r) << std::endl;
    std::cout << "prefactor: " << prefactor << std::endl;

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