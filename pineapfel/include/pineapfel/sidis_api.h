#pragma once

#include <apfel/apfelxx.h>

#include <map>
#include <string>
#include <vector>

namespace pineapfel {

// Mirror of apfel::SidisObjects that does not require including SIDIS.h.
// Having a single TU include SIDIS.h avoids ODR violations caused by
// APFEL++'s header-defined (non-inline) InitializeSIDIS functions.
struct SidisCoeffs {
    // LO
    apfel::DoubleObject<apfel::Operator>         C20qq;
    // NLO F2/FT
    apfel::DoubleObject<apfel::Operator>         C21qq;
    apfel::DoubleObject<apfel::Operator>         C21gq;
    apfel::DoubleObject<apfel::Operator>         C21qg;
    // NLO FL
    apfel::DoubleObject<apfel::Operator>         CL1qq;
    apfel::DoubleObject<apfel::Operator>         CL1gq;
    apfel::DoubleObject<apfel::Operator>         CL1qg;
    // NNLO FT non-singlet (nf-dependent)
    std::map<int, apfel::DoubleOperator>         C2Tns;
    // NNLO FT off-diagonal (nf-independent): "gq", "qg", "gg", "ps", "qbq",
    // "qpq1", "qpq2", "qpq3"
    std::map<std::string, apfel::DoubleOperator> C2T;
    // NNLO FL non-singlet (nf-dependent)
    std::map<int, apfel::DoubleOperator>         C2Lns;
    // NNLO FL off-diagonal (nf-independent)
    std::map<std::string, apfel::DoubleOperator> C2L;
};

// Thin wrapper around apfel::InitializeSIDIS (same x and z grid).
SidisCoeffs init_sidis(const apfel::Grid &g,
    const std::vector<double>            &thresholds);

// Mirror of apfel::SidisPolObjects (from SIDISpol.h).
// Same ODR-safe isolation as SidisCoeffs.
struct SidisPolCoeffs {
    // LO
    apfel::DoubleObject<apfel::Operator>         G10qq;
    // NLO G1
    apfel::DoubleObject<apfel::Operator>         G11qq;
    apfel::DoubleObject<apfel::Operator>         G11gq;
    apfel::DoubleObject<apfel::Operator>         G11qg;
    // NNLO G1 non-singlet (nf-dependent)
    std::map<int, apfel::DoubleOperator>         G12ns;
    // NNLO G1 gluon-in-quark (nf-dependent for polarized)
    std::map<int, apfel::DoubleOperator>         G12gq;
    // NNLO G1 quark-in-gluon (nf-dependent for polarized)
    std::map<int, apfel::DoubleOperator>         G12qg;
    // NNLO G1 nf-independent: "gg", "ps", "qbq", "qpq1", "qpq2", "qpq3"
    std::map<std::string, apfel::DoubleOperator> G12;
};

// Thin wrapper around apfel::InitializeSIDISpol (same x and z grid).
SidisPolCoeffs init_sidis_pol(const apfel::Grid &g,
    const std::vector<double>                   &thresholds);

} // namespace pineapfel
