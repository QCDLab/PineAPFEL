// This is the ONLY translation unit that includes <apfel/sidisbuilder.h>.

#include <pineapfel/sidis_api.h>

namespace pineapfel {

apfel::SidisNNLOObjects init_sidis_nnlo(const apfel::Grid &g,
    const std::vector<double>                             &thresholds,
    double                                                 int_eps) {
    return apfel::InitializeSidisObjects(g, g, thresholds, int_eps);
}

apfel::SidisNNLOObjects init_sidis_nnlo(const apfel::Grid &gx,
    const apfel::Grid                                     &gz,
    const std::vector<double>                             &thresholds,
    double                                                 int_eps) {
    return apfel::InitializeSidisObjects(gx, gz, thresholds, int_eps);
}

} // namespace pineapfel
