#include "catch.hpp"
#include "pseudopotential/sphericalHarmonics.h"
#include "qmfunctions/qmfunction_utils.h"



using namespace mrchem;

namespace density_tests {

TEST_CASE("sperical_harmonics", "[spherical_harmonics]") {
    const double prec = 1.0e-5;
    // const double thrs = 1.0e-12;

    SECTION("calc spherical harmonics") {
        double alpha = .5;
        OrbitalVector Phi;
        for (int l = 0; l < 5; l++) {
            for (int m = -l; m <= l; m++) {
                Phi.push_back(Orbital(SPIN::Paired));
            }
        }
        Phi.distribute();
        int iOrb = 0;
        for (int l = 0; l < 5; l++) {
            for (int m = -l; m <= l; m++) {
                auto fun = [l, m, alpha](const mrcpp::Coord<3> &r) {
                    double normr = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
                    switch (l) {
                        case 0:
                            return s0(r, normr) * std::exp(-alpha * normr);
                        case 1:
                            switch (m) {
                                case -1:
                                    return s1m1(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 1);
                                case 0:
                                    return s10(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 1);
                                case 1:
                                    return s11(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 1);
                            }
                        case 2:
                            switch (m) {
                                case -2:
                                    return s2m2(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 2);
                                case -1:
                                    return s2m1(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 2);
                                case 0:
                                    return s20(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 2);
                                case 1:
                                    return s21(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 2);
                                case 2:
                                    return s22(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 2);
                            }
                        case 3:
                            switch (m) {
                                case -3:
                                    return s3m3(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 3);
                                case -2:
                                    return s3m2(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 3);
                                case -1:
                                    return s3m1(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 3);
                                case 0:
                                    return s30(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 3);
                                case 1:                               
                                    return s31(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 3);
                                case 2:                               
                                    return s32(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 3);
                                case 3:                               
                                    return s33(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 3);
                            }
                        case 4:
                            switch (m) {
                                case -4:
                                    return s4m4(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 4);
                                case -3:
                                    return s4m3(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 4);
                                case -2:
                                    return s4m2(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 4);
                                case -1:
                                    return s4m1(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 4);
                                case 0:
                                    return s40(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 4);
                                case 1:
                                    return s41(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 4);
                                case 2:
                                    return s42(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 4);
                                case 3:
                                    return s43(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 4);
                                case 4:
                                    return s44(r, normr) * std::exp(-alpha * normr) / std::pow(normr, 4);
                            }
                    }
                    return 0.0;
                };
                if (mrcpp::mpi::my_orb(Phi[iOrb])) mrcpp::cplxfunc::project(Phi[iOrb], fun, NUMBER::Real, prec);
                iOrb++;
            }
        }
        Eigen::MatrixXd S(Phi.size(), Phi.size());
        for (int i = 0; i < Phi.size(); i++) {
            for (int j = 0; j < Phi.size(); j++) {
                // S(i, j) = qmfunction::dot(Phi[i], Phi[j]).real();
                S(i, j) = mrcpp::cplxfunc::dot(Phi[i], Phi[j]).real();
            }
        }
        Eigen::MatrixXd S_ref = Eigen::MatrixXd::Zero(Phi.size(), Phi.size());
        for (int i = 0; i < Phi.size(); i++) {
            S_ref(i, i) = 1.0 / (4.0 * alpha*alpha*alpha);
        }
        int iorb = 0;
        int jorb = 0;
        for (int l = 0; l < 5; l++) {
            for (int m = -l; m <= l; m++) {
                jorb = 0;
                for (int ll = 0; ll < 5; ll++) {
                    for (int mm = -ll; mm <= ll; mm++) {
                        if (l == ll && m == mm) {
                            std::cout << "l " << l << " m " << m << " ll " << ll << " mm " << mm << std::endl;
                            std::cout << iorb << " " << jorb << " " << S(iorb, jorb) << " " << S_ref(iorb, jorb) << std::endl;
                            REQUIRE(std::abs(S(iorb, jorb) - (S_ref(iorb, jorb))) < prec);
                        }
                        jorb++;
                    }
                }
                iorb++;
            }
        }
    }
}

} // namespace density_tests