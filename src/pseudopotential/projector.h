#pragma once
#include <mrchem.h>
#include "pseudopotential/sphericalHarmonics.h"

class ProjectorFunction: mrcpp::ComplexFunction {

public:
    ProjectorFunction(mrcpp::Coord<3> pos, double rl, int i, int l, int m, double prec);

    mrcpp::Coord<3> pos;
    int i;
    int l;
    int m;
    double rl;
    double prec;

private:
    /**
     * @brief Contains analytic form of projector.
    */
    double (*s)(const std::array<double, 3> &r, const double &normr);

    void switch_sperics(int l, int m);
    
};