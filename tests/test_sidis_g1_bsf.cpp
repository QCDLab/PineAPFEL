// test_sidis_g1_bsf.cpp
//
// End-to-end validation: PineAPFEL polarized SIDIS G₁ grid convolution
// vs APFEL++ BuildSidisG1 with toy PDFs/FFs.
//
// Setup:
//   - Toy PDFs/FFs (same function for both, Q-independent):
//       gluon:  xf(x) = x^{0.9}·(1-x)^5
//       quarks: xf(x) = x^{1.5}·(1-x)^3  (u,d,s,c,b and anti-quarks)
//   - Grid:    runcards/grid_sidis_nnlo_pol_g1.yaml  (12 bins, LO+NLO+NNLO)
//   - Theory:  runcards/sidis_g1_theory.yaml
//   - Ops:     runcards/sidis_g1_operator.yaml
//   - Order:   NNLO (PerturbativeOrder = 2, ch = -1)
//
// Convention notes (see fill.cpp comment for derivation):
//   PineAPFEL stores kernel weights scaled by (x_m·z_n)/(x_c·z_c·(4π)^n)
//   so that PineAPPL's convolution (which divides xf(x) by x) reproduces
//   exactly the APFEL++ BuildSidisG1 result G1_BSF(x_c, z_c, Q):
//
//     PineAPPL result = Σ_{Q²-nodes} G1_BSF(Q).Evaluate(x_c, z_c)
//
//   No additional coupling-factor or x/z rescaling is needed in this test.

#include <apfel/apfelxx.h>
#include <apfel/sidisbuilder.h>
#include <pineapfel/pineapfel.h>
#include <pineappl_capi.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <map>
#include <set>
#include <vector>

// ── Toy PDF/FF ───────────────────────────────────────────────────────────────
// Returns f(x) (NOT xf(x))
static double toy_f(int pid, double x) {
    if (pid == 21) return std::pow(x, -0.1) * std::pow(1.0 - x, 5.0);
    int apid = std::abs(pid);
    if (apid >= 1 && apid <= 5)
        return std::pow(x, 0.5) * std::pow(1.0 - x, 3.0);
    return 0.0;
}

// PineAPPL callback: returns x·f(x) (LHAPDF convention).
// PineAPPL then divides by x internally, yielding f(x).
// Same callback is used for both the PDF (x-slot) and FF (z-slot).
extern "C" double
xfx_callback(int32_t pid, double x, double /*q2*/, void * /*state*/) {
    return x * toy_f(pid, x);
}

// alpha_s callback backed by TabulateObject<double>
extern "C" double alphas_callback(double q2, void *state) {
    auto *as_tab = static_cast<apfel::TabulateObject<double> *>(state);
    return as_tab->Evaluate(std::sqrt(q2));
}

// ── Replicate fill.cpp's derive_q2_nodes ────────────────────────────────────
static std::vector<double> derive_q2_nodes(
    const std::vector<pineapfel::BinDef> &bins,
    const std::vector<double>            &thresholds,
    int                                   n_interm = 3) {
    std::set<double> q2_set;
    for (const auto &bin : bins) {
        double q2_lo = bin.lower[0], q2_hi = bin.upper[0];
        q2_set.insert(q2_lo);
        q2_set.insert(q2_hi);
        double log_lo = std::log(q2_lo), log_hi = std::log(q2_hi);
        for (int k = 1; k <= n_interm; k++) {
            double frac = static_cast<double>(k) / (n_interm + 1);
            q2_set.insert(std::exp(log_lo + frac * (log_hi - log_lo)));
        }
    }
    for (double thr : thresholds) {
        double q2t = thr * thr;
        if (!q2_set.empty() && q2t > *q2_set.begin() && q2t < *q2_set.rbegin())
            q2_set.insert(q2t);
    }
    return std::vector<double>(q2_set.begin(), q2_set.end());
}

int main() {
    std::printf("=== SIDIS G1 NNLO: PineAPFEL vs APFEL++ BuildSidisG1 ===\n\n");

    // ── Load runcards ────────────────────────────────────────────────────────
    auto grid_def =
        pineapfel::load_grid_def("runcards/grid_sidis_nnlo_pol_g1.yaml");
    auto theory = pineapfel::load_theory_card("runcards/sidis_g1_theory.yaml");
    auto op_card =
        pineapfel::load_operator_card("runcards/sidis_g1_operator.yaml");

    // ── alpha_s tabulation ───────────────────────────────────────────────────
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

    // ── APFEL++ grid from operator card ─────────────────────────────────────
    std::vector<apfel::SubGrid>   sgs;
    for (const auto &sg : op_card.xgrid)
        sgs.emplace_back(sg.n_knots, sg.x_min, sg.poly_degree);
    const apfel::Grid g{sgs};

    // ── Initialize SIDIS DoubleOperators (once, reused by BSF + grid build) ──
    std::printf(
        "Initializing SIDIS NNLO objects (this may take a moment)...\n");
    const auto sobj = pineapfel::init_sidis_nnlo(g,
        theory.quark_thresholds,
        op_card.sidis_int_eps);
    std::printf("Done.\n\n");

    // ── BuildSidisG1 reference with toy PDFs ────────────────────────────────
    // InPDFs / InFFs return xΔf(x) and zD(z) (LHAPDF convention).
    // Both PDF and FF use the same toy function.
    auto InDistFunc = [](double const &x,
                          double const & /*Q*/) -> std::map<int, double> {
        std::map<int, double> m;
        m[21] = x * toy_f(21, x);
        for (int pid = 1; pid <= 5; ++pid) {
            m[pid]  = x * toy_f(pid, x);
            m[-pid] = x * toy_f(-pid, x);
        }
        return m;
    };

    auto alphas_func = [&](double const &Q) -> double {
        return as_tab.Evaluate(Q);
    };

    // ch = -1 → full NNLO sum over all channels.
    auto G1_bsf = apfel::BuildSidisG1(sobj,
        InDistFunc,
        InDistFunc,
        alphas_func,
        theory.quark_thresholds,
        theory.pert_ord,
        -1);

    // ── Build PineAPFEL grid ─────────────────────────────────────────────────
    std::printf("Building PineAPFEL SIDIS G1 grid...\n");
    pineappl_grid *grid =
        pineapfel::build_grid(grid_def, theory, op_card, sobj);
    std::size_t nbins = pineappl_grid_bin_count(grid);
    std::printf("Grid built: %zu bins.\n\n", nbins);

    // ── PineAPPL convolution ─────────────────────────────────────────────────
    // Two convolutions: POL_PDF (x-slot) ⊗ UNPOL_FF (z-slot).
    // Both slots use the same toy callback.
    std::vector<double> predictions(nbins, 0.0);
    void               *pdfs_state[2] = {nullptr, nullptr};

    pineappl_grid_convolve(grid,
        xfx_callback,
        alphas_callback,
        pdfs_state,
        static_cast<void *>(&as_tab),
        nullptr, // all perturbative orders
        nullptr, // all channels
        nullptr, // no bin index filter
        0,       // nb_scales
        nullptr, // mu_scales
        predictions.data());

    // ── Q² nodes (must match fill.cpp logic) ────────────────────────────────
    auto q2_nodes = derive_q2_nodes(grid_def.bins, theory.quark_thresholds);
    std::printf("  Q² nodes: %zu\n\n", q2_nodes.size());

    // ── Comparison ───────────────────────────────────────────────────────────
    // PineAPFEL result == Σ_{Q²-nodes} G1_BSF(Q).Evaluate(x_c, z_c)
    //
    // Tolerance: 5%.  The small residual deviation comes from Lagrange
    // interpolation on the discrete grid; this matches test_f2nc_lhapdf.cpp.
    const double tolerance = 5e-2;

    std::printf("  %-8s  %-8s  %-14s  %-14s  %-10s  %s\n",
        "x_c",
        "z_c",
        "PineAPFEL",
        "BSF",
        "rel_diff",
        "");

    int failures = 0;
    for (std::size_t ibin = 0; ibin < nbins; ibin++) {
        double x_c = std::sqrt(
            grid_def.bins[ibin].lower[1] * grid_def.bins[ibin].upper[1]);
        double z_c = std::sqrt(
            grid_def.bins[ibin].lower[2] * grid_def.bins[ibin].upper[2]);

        // Sum BSF reference over the same Q² nodes that PineAPPL sums over.
        double ref = 0.0;
        for (double q2 : q2_nodes) {
            double Q  = std::sqrt(q2);
            ref      += G1_bsf(Q).Evaluate(x_c, z_c);
        }

        double pred    = predictions[ibin];
        double rel_err = (std::abs(ref) > 1e-30)
                             ? std::abs(pred - ref) / std::abs(ref)
                             : std::abs(pred);

        bool   ok      = (rel_err < tolerance);
        if (!ok) failures++;

        std::printf("  x=%.3e  z=%.3e  pred=%.6e  ref=%.6e  rel=%.2e  %s\n",
            x_c,
            z_c,
            pred,
            ref,
            rel_err,
            ok ? "OK" : "FAIL");
    }

    std::printf("\n=== Summary: %d / %zu bins FAILED ===\n\n", failures, nbins);

    pineappl_grid_delete(grid);
    return failures > 0 ? 1 : 0;
}
