#pragma once

#include <apfel/apfelxx.h>

#include <map>
#include <vector>

namespace pineapfel {

// Mirror of apfel::SidisObjects that does not require including SIDIS.h.
// Having a single TU include SIDIS.h avoids ODR violations caused by
// APFEL++'s header-defined (non-inline) InitializeSIDIS functions.
struct SidisCoeffs {
    apfel::DoubleObject<apfel::Operator>                C20qq;
    apfel::DoubleObject<apfel::Operator>                C21qq;
    apfel::DoubleObject<apfel::Operator>                C21gq;
    apfel::DoubleObject<apfel::Operator>                C21qg;
    apfel::DoubleObject<apfel::Operator>                CL1qq;
    apfel::DoubleObject<apfel::Operator>                CL1gq;
    apfel::DoubleObject<apfel::Operator>                CL1qg;
    std::map<int, apfel::DoubleObject<apfel::Operator>> C22qq;
};

// Thin wrapper around apfel::InitializeSIDIS (same x and z grid).
SidisCoeffs init_sidis(const apfel::Grid &g,
    const std::vector<double>            &thresholds);

} // namespace pineapfel
