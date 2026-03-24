// test_sidis_g1_bsf.cpp
//
// PineAPFEL polarised SIDIS G₁ (LO+NLO+NNLO) grid vs the same discrete
// DoubleOperator+z-quadrature reference used in test_grid_vs_apfelxx (sidis_helper).
//
// BuildSidisG1(...).Evaluate(x,z) summed over the PineAPPL Q² node list is *not*
// equal to PineAPPL's convolution of the imported subgrid: the grid stores
// per-node coefficient slices that PineAPPL sums with PDF/FF interpolation, while
// BuildSidisG1 assembles operators on fixed grids first (different discretisation).
// This test therefore uses compute_sidis_pol_reference + compute_sidis_g1_nnlo_reference.
//
// Fast setup: two bins, operator_test.yaml (same as TEST 6/8/12). NNLO object init
// is still ~O(N²) on the SIDIS grid (~2–3 minutes on a typical laptop).

#include <apfel/apfelxx.h>
#include <pineapfel/grid_gen.h>
#include <pineapfel/pineapfel.h>
#include <pineappl_capi.h>

#include "sidis_helper.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <set>
#include <vector>

static double toy_f(int pid, double x) {
    if (pid == 21) return std::pow(x, -0.1) * std::pow(1.0 - x, 5.0);
    int apid = std::abs(pid);
    if (apid >= 1 && apid <= 5)
        return std::pow(x, 0.5) * std::pow(1.0 - x, 3.0);
    return 0.0;
}

extern "C" double
xfx_callback(int32_t pid, double x, double /*q2*/, void * /*state*/) {
    return x * toy_f(pid, x);
}

extern "C" double alphas_callback(double q2, void *state) {
    auto *as_tab = static_cast<apfel::TabulateObject<double> *>(state);
    return as_tab->Evaluate(std::sqrt(q2));
}

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

static pineapfel::GridDef make_fast_sidis_g1_grid_def() {
    pineapfel::GridDef def;
    def.process            = pineapfel::ProcessType::SIDIS;
    def.observable         = pineapfel::Observable::F2;
    def.current            = pineapfel::Current::NC;
    def.cc_sign            = pineapfel::CCSign::Plus;
    def.mass_scheme        = pineapfel::MassScheme::ZM;
    def.pid_basis          = PINEAPPL_PID_BASIS_PDG;
    def.hadron_pids        = {2212, 211};
    def.convolution_types  = {PINEAPPL_CONV_TYPE_POL_PDF,
        PINEAPPL_CONV_TYPE_UNPOL_FF};
    def.orders             = {{0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0},
        {2, 0, 0, 0, 0}};
    def.bins               = {
        {{4.0, 0.008, 0.2}, {25.0, 0.04, 0.85}},
        {{25.0, 0.06, 0.2}, {200.0, 0.25, 0.85}},
    };
    def.normalizations     = {1.0, 1.0};
    return def;
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("=== SIDIS G1 NNLO: PineAPFEL vs discrete-operator reference ===\n\n");

    pineapfel::GridDef grid_def = make_fast_sidis_g1_grid_def();
    auto               theory =
        pineapfel::load_theory_card("runcards/sidis_g1_theory.yaml");
    auto op_card =
        pineapfel::load_operator_card("runcards/operator_test.yaml");

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

    std::vector<apfel::SubGrid> sgs;
    for (const auto &sg : op_card.xgrid)
        sgs.emplace_back(sg.n_knots, sg.x_min, sg.poly_degree);
    const apfel::Grid g{sgs};

    std::printf(
        "Initializing SIDIS NNLO objects (this may take a moment)...\n");
    const auto sobj = pineapfel::init_sidis_nnlo(g,
        theory.quark_thresholds,
        op_card.sidis_int_eps);
    std::printf("Done.\n\n");

    // ------------------------------------------------------------------
    // Sanity check: APFEL++ z integration uses Lagrange IntInterpolant.
    // APFEL++ does `dist_z.Integrate(z_lo,z_hi)` where dist_z is built
    // on the APFEL++ joint z-grid. We emulate the same integral by
    // summing IntInterpolant weights times dist_z evaluated at the
    // joint-grid nodes and require agreement.
    // ------------------------------------------------------------------
    {
        // Pick a kinematics point within the toy interpolation range.
        const double x_c  = 0.2;
        const double Q    = 10.0;
        const double z_lo = 0.2, z_hi = 0.85;

        const int nf = apfel::NF(Q, theory.quark_thresholds);

        // Toy inputs must provide all quark (±1..±nf) entries plus gluon.
        auto InpPDFsToy = [&](double const& x, double const& /*Q*/ )
                             -> std::map<int, double> {
            std::map<int, double> m;
            m[21] = x * toy_f(21, x);
            for (int j = -nf; j <= nf; j++) {
                if (j == 0) continue;
                m[j] = x * toy_f(j, x);
            }
            return m;
        };

        auto InFFsToy = [&](double const& z, double const& /*Q*/ )
                           -> std::map<int, double> {
            std::map<int, double> m;
            m[21] = z * toy_f(21, z);
            for (int j = -nf; j <= nf; j++) {
                if (j == 0) continue;
                m[j] = z * toy_f(j, z);
            }
            return m;
        };

        auto AsToy = [&](double const& mu) { return as_tab.Evaluate(mu); };

        // LO, polarized G1, channel=0 uses DoubleIdentity only.
        auto G1_lo = apfel::BuildSidisG1(sobj, InpPDFsToy, InFFsToy, AsToy,
            theory.quark_thresholds, 0, 0);

        const auto dist_z = G1_lo(Q).Evaluate1(x_c);
        const double apfel_int = dist_z.Integrate(z_lo, z_hi);

        const apfel::LagrangeInterpolator li{g};
        const auto &z_joint_nodes = g.GetJointGrid().GetGrid();
        const auto boundsa = li.SumBounds(z_lo, g.GetJointGrid());
        const auto boundsb = li.SumBounds(z_hi, g.GetJointGrid());

        double emu_int = 0.0;
        for (int beta = boundsa[0]; beta < boundsb[1]; beta++) {
            const double w = li.IntInterpolant(beta, z_lo, z_hi, g.GetJointGrid());
            const double z_node = z_joint_nodes[beta];
            emu_int += w * dist_z.Evaluate(z_node);
        }

        const double rel_err =
            (std::abs(apfel_int) > 0.0) ? std::abs(emu_int - apfel_int) / std::abs(apfel_int)
                                          : std::abs(emu_int - apfel_int);
        const double tol = 1e-10;
        if (rel_err > tol) {
            std::cerr << "IntInterpolant sanity check FAILED: apfel=" << apfel_int
                      << " emu=" << emu_int << " rel_err=" << rel_err << "\n";
            return 1;
        }
        std::printf("IntInterpolant sanity check: OK (rel_err=%.2e)\n\n", rel_err);
    }

    // ------------------------------------------------------------------
    // Focused diagnostic: LO point-mode with sidis_mode=bsf_exact should
    // reproduce APFEL++'s `BuildSidisG1(..., channel=-1, PerturbativeOrder=0)`
    // after `Evaluate1(x_c).Integrate(z_lo, z_hi)`.
    // ------------------------------------------------------------------
    {
        pineapfel::GridDef grid_def_point;
        grid_def_point.process            = pineapfel::ProcessType::SIDIS;
        grid_def_point.observable         = pineapfel::Observable::F2;
        grid_def_point.current            = pineapfel::Current::NC;
        grid_def_point.cc_sign            = pineapfel::CCSign::Plus;
        grid_def_point.mass_scheme        = pineapfel::MassScheme::ZM;
        grid_def_point.pid_basis          = PINEAPPL_PID_BASIS_PDG;
        grid_def_point.hadron_pids        = {2212, 211};
        grid_def_point.convolution_types  = {PINEAPPL_CONV_TYPE_POL_PDF,
            PINEAPPL_CONV_TYPE_UNPOL_FF};
        grid_def_point.orders = {{0, 0, 0, 0, 0}}; // LO only
        grid_def_point.bins   = {
            {{1.25001, 0.0052, 0.2}, {1.25001, 0.0052, 0.85}}};
        grid_def_point.normalizations = {1.0, 1.0};

        const double Q    = std::sqrt(1.25001);
        const double x_c  = 0.0052;
        const double z_lo = 0.2;
        const double z_hi = 0.85;
        const int    nf   = apfel::NF(Q, theory.quark_thresholds);

        auto InpPDFsToy = [&](double const& x, double const& /*Q*/) {
            std::map<int, double> m;
            m[21] = (x > 1.0) ? 0.0 : x * toy_f(21, x);
            for (int j = -nf; j <= nf; j++) {
                if (j == 0) continue;
                m[j] = (x > 1.0) ? 0.0 : x * toy_f(j, x);
            }
            return m;
        };
        auto InFFsToy = [&](double const& z, double const& /*Q*/) {
            std::map<int, double> m;
            m[21] = (z > 1.0) ? 0.0 : z * toy_f(21, z);
            for (int j = -nf; j <= nf; j++) {
                if (j == 0) continue;
                m[j] = (z > 1.0) ? 0.0 : z * toy_f(j, z);
            }
            return m;
        };

        auto AsToy = [&](double const& mu) { return as_tab.Evaluate(mu); };

        auto G1_lo_builder = apfel::BuildSidisG1(sobj,
            InpPDFsToy,
            InFFsToy,
            AsToy,
            theory.quark_thresholds,
            /*PerturbativeOrder=*/0,
            /*channel=*/-1);

        const double apfel_int = G1_lo_builder(Q).Evaluate1(x_c).Integrate(z_lo, z_hi);

        const char *modes[] = {"legacy", "bsf_exact"};
        const double tol = 1e-10;
        int any_fail = 0;
        for (const char *mode : modes) {
            auto op_card_mode = op_card;
            op_card_mode.sidis_mode = mode;

            pineappl_grid *grid_point =
                pineapfel::build_grid(grid_def_point,
                    theory,
                    op_card_mode,
                    sobj);

            std::vector<double> pred_point(
                pineappl_grid_bin_count(grid_point), 0.0);
            void *pdfs_state[2] = {nullptr, nullptr};

            pineappl_grid_convolve(grid_point,
                xfx_callback,
                alphas_callback,
                pdfs_state,
                static_cast<void *>(&as_tab),
                nullptr,
                nullptr,
                nullptr,
                0,
                nullptr,
                pred_point.data());

            pineappl_grid_delete(grid_point);

            const double rel =
                (std::abs(apfel_int) > 0.0)
                    ? std::abs(pred_point[0] - apfel_int) / std::abs(apfel_int)
                    : std::abs(pred_point[0] - apfel_int);
            if (rel > tol) {
                any_fail = 1;
                std::cerr
                    << "LO point-mode diagnostic FAILED for mode=" << mode
                    << " pine=" << pred_point[0] << " apfel=" << apfel_int
                    << " rel=" << rel << "\n";
            } else {
                std::printf("LO point-mode diagnostic: mode=%s OK (rel_err=%.2e)\n\n",
                    mode,
                    rel);
            }
        }
        if (any_fail) return 1;
    }

    std::printf("Building PineAPFEL SIDIS G1 grid...\n");
    pineappl_grid *grid =
        pineapfel::build_grid(grid_def, theory, op_card, sobj);
    std::size_t    nbins = pineappl_grid_bin_count(grid);
    std::printf("Grid built: %zu bins.\n\n", nbins);

    std::vector<double> predictions(nbins, 0.0);
    void               *pdfs_state[2] = {nullptr, nullptr};

    pineappl_grid_convolve(grid,
        xfx_callback,
        alphas_callback,
        pdfs_state,
        static_cast<void *>(&as_tab),
        nullptr,
        nullptr,
        nullptr,
        0,
        nullptr,
        predictions.data());

    auto                        q2_nodes = derive_q2_nodes(
        grid_def.bins, theory.quark_thresholds);
    std::vector<std::vector<double>> bin_x, bin_z;
    bin_x.reserve(nbins);
    bin_z.reserve(nbins);
    double q2_top = 0.0;
    for (const auto &b : grid_def.bins) {
        bin_x.push_back({b.lower[1], b.upper[1]});
        bin_z.push_back({b.lower[2], b.upper[2]});
        q2_top = std::max(q2_top, b.upper[0]);
    }
    const int nf_max = apfel::NF(std::sqrt(q2_top), theory.quark_thresholds);

    auto toy_wrap = [](int p, double x) { return toy_f(p, x); };
    auto alphas_wrap =
        [&](double Q) { return as_tab.Evaluate(Q); };

    auto ref_lo_nlo = compute_sidis_pol_reference(g,
        sobj,
        theory.quark_thresholds,
        q2_nodes,
        bin_x,
        bin_z,
        1,
        toy_wrap,
        alphas_wrap);
    auto ref_nnlo = compute_sidis_g1_nnlo_reference(g,
        sobj,
        theory.quark_thresholds,
        nf_max,
        q2_nodes,
        bin_x,
        bin_z,
        toy_wrap,
        alphas_wrap);

    const double tolerance = 1e-5;

    std::printf("  Q² nodes: %zu  nf_max=%d\n\n", q2_nodes.size(), nf_max);
    std::printf("  %-8s  %-14s  %-14s  %-14s  %-10s  %s\n",
        "x_c",
        "z_bin",
        "PineAPFEL",
        "ref",
        "rel_diff",
        "");

    int failures = 0;
    for (std::size_t ibin = 0; ibin < nbins; ibin++) {
        double x_c  = std::sqrt(
            grid_def.bins[ibin].lower[1] * grid_def.bins[ibin].upper[1]);
        double z_lo = grid_def.bins[ibin].lower[2];
        double z_hi = grid_def.bins[ibin].upper[2];

        double ref =
            ref_lo_nlo[ibin] + ref_nnlo[ibin];
        double pred    = predictions[ibin];
        double rel_err = (std::abs(ref) > 1e-30)
                             ? std::abs(pred - ref) / std::abs(ref)
                             : std::abs(pred);

        bool ok = (rel_err < tolerance);
        if (!ok) failures++;

        std::printf(
            "  x=%.3e  z=[%.2f,%.2f]  pred=%.6e  ref=%.6e  rel=%.2e  %s\n",
            x_c,
            z_lo,
            z_hi,
            pred,
            ref,
            rel_err,
            ok ? "OK" : "FAIL");
    }

    std::printf("\n=== Summary: %d / %zu bins FAILED ===\n\n", failures, nbins);

    pineappl_grid_delete(grid);
    return failures > 0 ? 1 : 0;
}
