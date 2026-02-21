// This is the ONLY translation unit that includes <apfel/SIDISpol.h>
// (which in turn includes <apfel/SIDIS.h>).  Both headers define non-inline
// functions, so including them in multiple TUs causes ODR violations at link
// time.  All other files use pineapfel::init_sidis() / init_sidis_pol().

#include <sidis_api.h>

// SIDISpol.h already includes SIDIS.h — include only the outermost header.
#include <apfel/SIDISpol.h>

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

SidisPolCoeffs init_sidis_pol(const apfel::Grid &g,
    const std::vector<double>                   &thresholds) {
    auto           s = apfel::InitializeSIDISpol(g, thresholds);
    SidisPolCoeffs c;
    c.G10qq = s.G10qq;
    c.G11qq = s.G11qq;
    c.G11gq = s.G11gq;
    c.G11qg = s.G11qg;
    c.G12qq = s.G12qq;
    return c;
}

} // namespace pineapfel
