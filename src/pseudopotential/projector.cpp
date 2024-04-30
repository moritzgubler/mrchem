
#include "pseudopotential/projector.h"
#include <math.h>

Projector::Projector(Vector3d pos, double rl, int i, int l, int m) {
    this->pos = pos;
    this->rl = rl;
    this->i = i;
    this->l = l;
    this->m = m;
    // switch_sperics(l, m);
    double prefactor = std::sqrt(2.0) / (std::pow(rl, l + (4 * i - 2) / 2) * std::sqrt(tgamma( l + (4 * i - 2) / 2 )) );

    auto project_analytic = [this, prefactor](const Eigen::Vector3d &r) {
        double normr = r.norm();
        return prefactor * std::pow(normr, this->l + 2 * (this->i - 1)) * this->s(r, normr);
    };
    auto op = (*this);
    

}

void Projector::switch_sperics(int l, int m){
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