// test_f2nc_lhapdf.cpp
//
// End-to-end validation: PineAPFEL F2 NC unpolarized (grid convolution)
// vs APFEL++ BuildStructureFunctions with toy PDFs.
//
// Setup:
//   - Toy PDFs : g(x) = x^{-0.1}(1-x)^5, q(x) = x^{0.5}(1-x)^3 for u,d,s,c,b
//   - x-grid   : {SubGrid{100,1e-5,3}, SubGrid{60,1e-1,3},
//                 SubGrid{50,6e-1,3},  SubGrid{50,8e-1,3}}
//   - Order    : NNLO (PerturbativeOrder = 2)

#include <apfel/apfelxx.h>
#include <pineapfel/fill.h>
#include <pineappl_capi.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <map>
#include <set>
#include <vector>

// ── Toy PDFs ────────────────────────────────────────────────────────────────
// toy_f returns f(x), NOT x*f(x)
static double toy_f(int pid, double x) {
    if (pid == 21) return std::pow(x, -0.1) * std::pow(1.0 - x, 5.0);
    int apid = std::abs(pid);
    if (apid >= 1 && apid <= 5)
        return std::pow(x, 0.5) * std::pow(1.0 - x, 3.0);
    return 0.0;
}

// PineAPPL callback: returns x*f(x) (LHAPDF convention)
extern "C" double
xfx_callback(int32_t pid, double x, double /*q2*/, void * /*state*/) {
    return x * toy_f(pid, x);
}

// alpha_s callback backed by TabulateObject<double>
extern "C" double alphas_callback(double q2, void *state) {
    auto *as_tab = static_cast<apfel::TabulateObject<double> *>(state);
    return as_tab->Evaluate(std::sqrt(q2));
}

int main() {
    std::printf("=== F2 NC Unpol: PineAPFEL vs APFEL++ BSF (toy PDFs) ===\n\n");

    // ── Theory setup ────────────────────────────────────────────────────────
    auto grid_def = pineapfel::load_grid_def("runcards/grid_dis_unp_nc.yaml");
    auto theory   = pineapfel::load_theory_card("runcards/theory.yaml");
    auto op_card =
        pineapfel::load_operator_card("runcards/operator_apfelxx.yaml");

    // ── alpha_s tabulation ──────────────────────────────────────────────────
    apfel::AlphaQCD               a{theory.alpha_qcd_ref,
        theory.q_ref,
        theory.quark_thresholds,
        theory.pert_ord};
    const auto                   &tabp = op_card.tabulation;
    apfel::TabulateObject<double> as_tab{a,
        tabp.n_points,
        tabp.q_min,
        static_cast<double>(tabp.n_steps),
        tabp.interp_degree};

    // ── APFEL++ BSF reference ───────────────────────────────────────────────
    // x-grid matching operator_apfelxx.yaml
    std::vector<apfel::SubGrid>   sgs;
    for (const auto &sg : op_card.xgrid)
        sgs.emplace_back(sg.n_knots, sg.x_min, sg.poly_degree);
    const apfel::Grid g{sgs};

    // Evolution-basis InDistFunc (xf convention, matching test 4 in
    // test_grid_vs_apfelxx.cpp).  For 5 equal-mass light quarks:
    //   GLUON(0)  = x * g(x)
    //   SIGMA(1)  = 2*5 * x*f(x) = 10 * x*f(x)
    //   T35(11)   = 10 * x*f(x)   (5 flavours, no top)
    auto              dist_func =
        [](int const &i, double const &x, double const & /*Q*/) -> double {
        double f     = std::pow(x, 0.5) * std::pow(1.0 - x, 3.0);
        double g_val = std::pow(x, -0.1) * std::pow(1.0 - x, 5.0);
        if (i == 0) return x * g_val;     // GLUON (xf convention)
        if (i == 1) return 10.0 * x * f;  // SIGMA (xf convention)
        if (i == 11) return 10.0 * x * f; // T35   (xf convention)
        return 0.0;
    };

    auto alphas_func = [&](double const &Q) -> double {
        return as_tab.Evaluate(Q);
    };
    auto couplings_func = [&](double const &Q) -> std::vector<double> {
        return apfel::ElectroWeakCharges(Q, false);
    };

    auto F2Obj = apfel::InitializeF2NCObjectsZM(g, theory.quark_thresholds);
    auto F2    = apfel::BuildStructureFunctions(F2Obj,
        dist_func,
        theory.pert_ord,
        alphas_func,
        couplings_func);

    // ── PineAPFEL grid ──────────────────────────────────────────────────────
    pineappl_grid      *grid = pineapfel::build_grid(grid_def, theory, op_card);
    std::size_t         nbins = pineappl_grid_bin_count(grid);

    std::vector<double> predictions(nbins, 0.0);
    pineappl_grid_convolve_with_one(grid,
        2212,
        xfx_callback,
        alphas_callback,
        static_cast<void *>(&as_tab),
        nullptr, // all perturbative orders
        nullptr, // all channels
        1.0,     // xi_ren
        1.0,     // xi_fac
        predictions.data());

    // ── Derive Q² nodes (same logic as fill.cpp derive_q2_nodes) ───────────
    const int        n_intermediate = 3;
    std::set<double> q2_set;
    for (const auto &bin : grid_def.bins) {
        double q2_lo = bin.lower[0], q2_hi = bin.upper[0];
        q2_set.insert(q2_lo);
        q2_set.insert(q2_hi);
        double log_lo = std::log(q2_lo), log_hi = std::log(q2_hi);
        for (int k = 1; k <= n_intermediate; k++) {
            double frac = static_cast<double>(k) / (n_intermediate + 1);
            q2_set.insert(std::exp(log_lo + frac * (log_hi - log_lo)));
        }
    }
    for (double thr : theory.quark_thresholds) {
        double q2t = thr * thr;
        if (!q2_set.empty() && q2t > *q2_set.begin() && q2t < *q2_set.rbegin())
            q2_set.insert(q2t);
    }
    std::vector<double> q2_nodes(q2_set.begin(), q2_set.end());

    std::printf("  Q² nodes used: %zu\n", q2_nodes.size());
    for (double q2 : q2_nodes)
        std::printf("    Q²=%.4e  Q=%.4f\n", q2, std::sqrt(q2));
    std::printf("\n");

    // ── Comparison ──────────────────────────────────────────────────────────
    // 5% tolerance: PineAPPL uses interpolated subgrids while BSF computes
    // exact convolutions.
    const double tolerance = 5e-2;

    std::printf("  %-10s  %-14s  %-14s  %-10s  %s\n",
        "x",
        "PineAPFEL",
        "BSF",
        "rel_diff",
        "");

    int failures = 0;
    for (std::size_t ibin = 0; ibin < nbins; ibin++) {
        double x_lo = grid_def.bins[ibin].lower.back();
        double x_hi = grid_def.bins[ibin].upper.back();
        double x_c  = std::sqrt(x_lo * x_hi);

        // Sum BSF reference over the same q2_nodes that PineAPPL sums over.
        double ref  = 0.0;
        for (double q2 : q2_nodes) ref += F2.at(0).Evaluate(x_c, std::sqrt(q2));

        double pred = predictions[ibin];
        double rel  = (std::abs(ref) > 1e-30)
                          ? std::abs(pred - ref) / std::abs(ref)
                          : std::abs(pred);

        bool   ok   = (rel < tolerance);
        if (!ok) failures++;

        std::printf("  x=%.2e  pred=%.6e  ref=%.6e  rel_diff=%.2e  %s\n",
            x_c,
            pred,
            ref,
            rel,
            ok ? "OK" : "FAIL");
    }

    std::printf("\n=== Summary: %d / %zu bins FAILED ===\n\n", failures, nbins);

    pineappl_grid_delete(grid);

    return failures > 0 ? 1 : 0;
}
