// This is the ONLY translation unit that includes <apfel/SIDIS.h>.
// InitializeSIDIS is defined (not declared) there, so including it in
// multiple TUs causes ODR violations at link time.  All other files that
// need SIDIS coefficient functions use pineapfel::init_sidis() instead.

#include <sidis_api.h>

#include <apfel/SIDIS.h>

namespace pineapfel {

SidisCoeffs init_sidis(const apfel::Grid &g,
    const std::vector<double>            &thresholds) {
    auto        s = apfel::InitializeSIDIS(g, thresholds);
    SidisCoeffs c;
    c.C20qq = s.C20qq;
    c.C21qq = s.C21qq;
    c.C21gq = s.C21gq;
    c.C21qg = s.C21qg;
    c.CL1qq = s.CL1qq;
    c.CL1gq = s.CL1gq;
    c.CL1qg = s.CL1qg;
    c.C22qq = s.C22qq;
    return c;
}

} // namespace pineapfel
