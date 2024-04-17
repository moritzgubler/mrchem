#ifndef SURFACEFORCE_H
#define SURFACEFORCE_H

#include <vector>
#include "qmfunctions/Orbital.h"
#include "chemistry/Molecule.h"

namespace surface_force {

// Function declaration
std::vector<double> surface_forces(mrchem::Molecule &mol, mrchem::OrbitalVector &Phi, double prec);

} // namespace surface_force

#endif // SURFACEFORCE_H

