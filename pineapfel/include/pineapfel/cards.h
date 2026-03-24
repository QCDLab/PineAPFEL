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
    // Optional z-grid for SIDIS. If omitted, xgrid is reused for z.
    std::vector<SubGridDef> zgrid;
    TabulationParams        tabulation;
    std::vector<double>     xi;
    // Relative integration precision for SIDIS NNLO coefficient functions.
    // Controls the cost of InitializeSidisObjects.
    // Default 1e-3 gives accurate results; tests may 1e-1 for ~100× speedup
    // while keeping PineAPPL vs APFEL++ exact (both see the same operator).
    double                  sidis_int_eps = 1e-3;
    // Controls the SIDIS exact (BSF) grid filling implementation.
    // - "legacy": use the previous external z quadrature approach.
    // - "bsf_exact": use APFEL++-consistent z integration via IntInterpolant.
    std::string             sidis_mode = "legacy";
    // SIDIS fill controls for Q^2/z integration strategy.
    // q2_n_intermediate:
    //   number of geometric interior Q^2 nodes per bin edge interval.
    // q2_include_thresholds:
    //   insert heavy-quark threshold Q^2 points in global node set.
    // q2_use_bin_centers_only:
    //   ignore edge/intermediate strategy and use only bin geometric centres.
    // z_quad_subdivisions:
    //   split each z-bin into this many sub-intervals, each integrated
    //   with the internal 8-point Gauss-Legendre rule.
    int                     sidis_q2_n_intermediate      = 3;
    bool                    sidis_q2_include_thresholds  = true;
    bool                    sidis_q2_use_bin_centers_only = false;
    int                     sidis_z_quad_subdivisions    = 1;
};

TheoryCard   load_theory_card(const std::string &path);
OperatorCard load_operator_card(const std::string &path);

} // namespace pineapfel
