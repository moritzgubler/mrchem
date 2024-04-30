#include <mrchem.h>
#include "pseudopotential/sphericalHarmonics.h"

class Projector: mrcpp::ComplexFunction {

public:
    Projector(Vector3d pos, double rl, int i, int l, int m);

    Vector3d pos;
    int i;
    int l;
    int m;
    double rl;

private: 
    double (*s)(const Eigen::Vector3d &r, const double &normr);

    void switch_sperics(int l, int m);
    
};