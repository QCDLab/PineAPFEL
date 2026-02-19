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
    bool                qed;
    double              alpha_qed_ref;
    std::vector<double> lepton_thresholds;
};

struct OperatorCard {
    std::vector<SubGridDef> xgrid;
    TabulationParams        tabulation;
    std::vector<double>     xi;
};

TheoryCard   load_theory_card(const std::string &path);
OperatorCard load_operator_card(const std::string &path);

} // namespace pineapfel
