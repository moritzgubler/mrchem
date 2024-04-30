
double s0(const Eigen::Vector3d &r, const double &normr) {
    return 0.5 * std::sqrt(1.0 / M_PI);
}

double s1m1(const Eigen::Vector3d &r, const double &normr) {
    return std::sqrt(3.0 / (4.0 * M_PI)) * r(1) / normr;
}

double s10(const Eigen::Vector3d &r, const double &normr) {
    return std::sqrt(3.0 / (4.0 * M_PI)) * r(2) / normr;
}

double s11(const Eigen::Vector3d &r, const double &normr) {
    return std::sqrt(3.0 / (4.0 * M_PI)) * r(0) / normr;
}

double s2m2(const Eigen::Vector3d &r, const double &normr) {
    return 0.5 * std::sqrt(15.0 / ( M_PI)) * r(0) * r(1) / ( normr * normr);
}

double s2m1(const Eigen::Vector3d &r, const double &normr) {
    return 0.5 * std::sqrt(15.0 / ( M_PI)) * r(1) * r(2) / (normr * normr);
}

double s20(const Eigen::Vector3d &r, const double &normr) {
    return 0.25 * std::sqrt(5.0 / ( M_PI)) * (3.0 * r(2) * r(2) - normr * normr) / (normr * normr);
}

double s21(const Eigen::Vector3d &r, const double &normr) {
    return 0.5 * std::sqrt(15.0 / ( M_PI)) * r(0) * r(2) / (normr * normr);
}

double s22(const Eigen::Vector3d &r, const double &normr) {
    return 0.25 * std::sqrt(15.0 / ( M_PI)) * (r(0) * r(0) - r(1) * r(1)) / (normr * normr);
}

double s3m3(const Eigen::Vector3d &r, const double &normr) {
    return 0.25 * std::sqrt(35.0 / ( 2 * M_PI)) * (r(1) * (3 * r(0) * r(0) - r(1) * r(1)) / (normr * normr * normr) );
}

double s3m2(const Eigen::Vector3d &r, const double &normr) {
    return 0.5 * std::sqrt(105.0 / ( M_PI)) * r(0) * r(1) * r(2) / (normr * normr * normr);
}

double s3m1(const Eigen::Vector3d &r, const double &normr) {
    return 0.25 * std::sqrt(21.0 / ( 2 * M_PI)) * r(1) * (5 * r(2) * r(2) - normr * normr) / (normr * normr * normr);
}

double s30(const Eigen::Vector3d &r, const double &normr) {
    return 0.25 * std::sqrt(7.0 / ( M_PI)) * (5 * r(2) * r(2) * r(3) - 3 * r(2) * normr * normr ) / (normr * normr * normr);
}

double s31(const Eigen::Vector3d &r, const double &normr) {
    return 0.25 * std::sqrt(21.0 / ( 2 * M_PI)) * r(0) * (5 * r(2) * r(2) - normr * normr) / (normr * normr * normr);
}

double s32(const Eigen::Vector3d &r, const double &normr) {
    return 0.25 * std::sqrt(105.0 / ( M_PI)) * (r(0) * r(0) - r(1) * r(1)) * r(2) / (normr * normr * normr);
}

double s33(const Eigen::Vector3d &r, const double &normr) {
    return 0.25 * std::sqrt(35.0 / ( 2 * M_PI)) * (r(0) * r(0) * r(0) - 3 * r(1) * r(1) * r(0)) / (normr * normr * normr);
}