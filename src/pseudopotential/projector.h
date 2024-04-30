#include <mrchem.h>
#include "pseudopotential/sphericalHarmonics.h"

class Projector: mrcpp::ComplexFunction {

public:
    Projector(Vector3d pos, double rl, int i, int l, int m, double prec);

    Vector3d pos;
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