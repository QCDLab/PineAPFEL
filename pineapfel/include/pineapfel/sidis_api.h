#pragma once

#include <apfel/apfelxx.h>
#include <apfel/sidisbuilder.h>

#include <vector>

namespace pineapfel {

// Returns pre-computed SIDIS DoubleOperator objects (LO + NLO + exact NNLO)
// for unpolarised FT, FL and polarised G1 structure functions.
// Both PDF and FF use the same grid g.
apfel::SidisNNLOObjects init_sidis_nnlo(const apfel::Grid &g,
    const std::vector<double>                             &thresholds,
    double                                                 int_eps = 1e-3);

// Overload with separate x and z grids (Zurich/APFEL++ benchmark setup).
apfel::SidisNNLOObjects init_sidis_nnlo(const apfel::Grid &gx,
    const apfel::Grid                                        &gz,
    const std::vector<double>                             &thresholds,
    double                                                 int_eps = 1e-3);

} // namespace pineapfel
