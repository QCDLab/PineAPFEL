// This is the ONLY translation unit that includes <apfel/sidisbuilder.h>.
// sidisbuilder.h uses #pragma once and only declares (no header-defined
// function bodies), so there are no ODR violations.  All other files access
// SIDIS coefficient functions through pineapfel::init_sidis_nnlo().

#include <pineapfel/sidis_api.h>

namespace pineapfel {

apfel::SidisNNLOObjects init_sidis_nnlo(const apfel::Grid &g,
    const std::vector<double>                             &thresholds,
    double                                                 int_eps) {
    return apfel::InitializeSidisObjects(g, g, thresholds, int_eps);
}

} // namespace pineapfel
