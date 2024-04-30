/**
 * @file sphericalHarmonics.h
 * @brief Contains the relevant, real spherical harmonics functions for the non local pseudopotential.
 * The naming of the functions is as follows: slm where l is the angular momentum and m is the magnetic quantum number. m stands for minus.
 * An example of the naming is s2m1 which is the spherical harmonic function with l=2 and m=-1.
 */

double s0(const std::array<double, 3> &r, const double &normr) {
    return 0.5 * std::sqrt(1.0 / M_PI);
}

double s1m1(const std::array<double, 3> &r, const double &normr) {
    return std::sqrt(3.0 / (4.0 * M_PI)) * r[1] / normr;
}

double s10(const std::array<double, 3> &r, const double &normr) {
    return std::sqrt(3.0 / (4.0 * M_PI)) * r[2] / normr;
}

double s11(const std::array<double, 3> &r, const double &normr) {
    return std::sqrt(3.0 / (4.0 * M_PI)) * r[0] / normr;
}

double s2m2(const std::array<double, 3> &r, const double &normr) {
    return 0.5 * std::sqrt(15.0 / ( M_PI)) * r[0] * r[1] / ( normr * normr);
}

double s2m1(const std::array<double, 3> &r, const double &normr) {
    return 0.5 * std::sqrt(15.0 / ( M_PI)) * r[1] * r[2] / (normr * normr);
}

double s20(const std::array<double, 3> &r, const double &normr) {
    return 0.25 * std::sqrt(5.0 / ( M_PI)) * (3.0 * r[2] * r[2] - normr * normr) / (normr * normr);
}

double s21(const std::array<double, 3> &r, const double &normr) {
    return 0.5 * std::sqrt(15.0 / ( M_PI)) * r[0] * r[2] / (normr * normr);
}

double s22(const std::array<double, 3> &r, const double &normr) {
    return 0.25 * std::sqrt(15.0 / ( M_PI)) * (r[0] * r[0] - r[1] * r[1]) / (normr * normr);
}

double s3m3(const std::array<double, 3> &r, const double &normr) {
    return 0.25 * std::sqrt(35.0 / ( 2 * M_PI)) * (r[1] * (3 * r[0] * r[0] - r[1] * r[1]) / (normr * normr * normr) );
}

double s3m2(const std::array<double, 3> &r, const double &normr) {
    return 0.5 * std::sqrt(105.0 / ( M_PI)) * r[0] * r[1] * r[2] / (normr * normr * normr);
}

double s3m1(const std::array<double, 3> &r, const double &normr) {
    return 0.25 * std::sqrt(21.0 / ( 2 * M_PI)) * r[1] * (5 * r[2] * r[2] - normr * normr) / (normr * normr * normr);
}

double s30(const std::array<double, 3> &r, const double &normr) {
    return 0.25 * std::sqrt(7.0 / ( M_PI)) * (5 * r[2] * r[2] * r[2] - 3 * r[2] * normr * normr ) / (normr * normr * normr);
}

double s31(const std::array<double, 3> &r, const double &normr) {
    return 0.25 * std::sqrt(21.0 / ( 2 * M_PI)) * r[0] * (5 * r[2] * r[2] - normr * normr) / (normr * normr * normr);
}

double s32(const std::array<double, 3> &r, const double &normr) {
    return 0.25 * std::sqrt(105.0 / ( M_PI)) * (r[0] * r[0] - r[1] * r[1]) * r[2] / (normr * normr * normr);
}

double s33(const std::array<double, 3> &r, const double &normr) {
    return 0.25 * std::sqrt(35.0 / ( 2 * M_PI)) * (r[0] * r[0] * r[0] - 3 * r[1] * r[1] * r[0]) / (normr * normr * normr);
}