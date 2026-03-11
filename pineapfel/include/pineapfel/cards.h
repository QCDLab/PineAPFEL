#pragma once

#include <string>
#include <vector>

namespace pineapfel {

struct SubGridDef {
    int    n_knots;
    double x_min;
    int    poly_degree;
};

struct TabulationParams {
    int    n_points;
    double q_min;
    int    n_steps;
    int    interp_degree;
};

struct TheoryCard {
    double              mu0;
    int                 pert_ord;
    double              q_ref;
    double              alpha_qcd_ref;
    std::vector<double> quark_thresholds;
    std::vector<int>    flavors;
    std::vector<double> ckm;
    bool                qed;
    double              alpha_qed_ref;
    std::vector<double> lepton_thresholds;
    // Heavy quark masses for massive scheme coefficient functions
    std::vector<double> heavy_quark_masses;
    // Tabulation parameters for massive scheme (APFEL++ defaults)
    int                 mass_nxi    = 150;
    double              mass_ximin  = 0.05;
    double              mass_ximax  = 10000.0;
    int                 mass_intdeg = 3;
    double              mass_lambda = 0.0005;
    int                 mass_imod   = 0;
};

struct OperatorCard {
    std::vector<SubGridDef> xgrid;
    TabulationParams        tabulation;
    std::vector<double>     xi;
};

TheoryCard   load_theory_card(const std::string &path);
OperatorCard load_operator_card(const std::string &path);

} // namespace pineapfel
