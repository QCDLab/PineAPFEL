#include "sidis_helper.h"
#include <apfel/apfelxx.h>
#include <pineapfel.h>
#include <pineappl_capi.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <vector>

// ----------------------------------------------------------------
// Toy PDF: returns f(x), NOT x*f(x)
// ----------------------------------------------------------------
static double toy_f(int pid, double x) {
    if (pid == 21) return std::pow(x, -0.1) * std::pow(1.0 - x, 5.0);
    int apid = std::abs(pid);
    if (apid >= 1 && apid <= 5)
        return std::pow(x, 0.5) * std::pow(1.0 - x, 3.0);
    return 0.0;
}

// PineAPPL callback: returns x*f(x)
extern "C" double
xfx_callback(int32_t pid, double x, double /*q2*/, void * /*state*/) {
    return x * toy_f(pid, x);
}

// Alpha_s callback
extern "C" double alphas_callback(double q2, void *state) {
    auto *as_tab = static_cast<apfel::TabulateObject<double> *>(state);
    return as_tab->Evaluate(std::sqrt(q2));
}

// ----------------------------------------------------------------
// Reproduce build_channel_operator from fill.cpp
// ----------------------------------------------------------------
static apfel::Operator build_channel_operator(
    const pineapfel::ChannelDef          &channel,
    const std::map<int, apfel::Operator> &ops,
    const std::vector<double>            &charges,
    int                                   nf) {
    const apfel::Operator &CNS    = ops.at(apfel::DISNCBasis::CNS);
    const apfel::Operator &CS     = ops.at(apfel::DISNCBasis::CS);
    const apfel::Operator &CG     = ops.at(apfel::DISNCBasis::CG);

    double                 sum_ch = 0;
    for (int i = 0; i < nf && i < static_cast<int>(charges.size()); i++)
        sum_ch += charges[i];

    bool is_gluon = false;
    for (const auto &combo : channel.pid_combinations) {
        if (combo.size() == 1 && combo[0] == 21) {
            is_gluon = true;
            break;
        }
    }

    if (is_gluon) { return sum_ch * CG; }

    int quark_pid = 0;
    for (const auto &combo : channel.pid_combinations) {
        if (combo.size() == 1 && combo[0] > 0) {
            quark_pid = combo[0];
            break;
        }
    }
    if (quark_pid == 0) {
        for (const auto &combo : channel.pid_combinations) {
            if (combo.size() == 1 && combo[0] < 0) {
                quark_pid = -combo[0];
                break;
            }
        }
    }

    int    q_idx  = quark_pid - 1;
    double e_q_sq = (q_idx >= 0 && q_idx < static_cast<int>(charges.size()))
                        ? charges[q_idx]
                        : 0.0;

    return e_q_sq * CNS + (sum_ch / 6.0) * (CS - CNS);
}

// ----------------------------------------------------------------
// Reproduce derive_q2_nodes from fill.cpp
// ----------------------------------------------------------------
static std::vector<double> derive_q2_nodes(
    const std::vector<pineapfel::BinDef> &bins,
    const std::vector<double>            &thresholds,
    int                                   n_intermediate = 3) {
    std::set<double> q2_set;
    for (const auto &bin : bins) {
        double q2_lo = bin.lower[0], q2_hi = bin.upper[0];
        q2_set.insert(q2_lo);
        q2_set.insert(q2_hi);
        if (n_intermediate > 0) {
            double log_lo = std::log(q2_lo), log_hi = std::log(q2_hi);
            for (int i = 1; i <= n_intermediate; i++) {
                double frac = static_cast<double>(i) / (n_intermediate + 1);
                q2_set.insert(std::exp(log_lo + frac * (log_hi - log_lo)));
            }
        }
    }
    double q2_min = *q2_set.begin(), q2_max = *q2_set.rbegin();
    for (double thr : thresholds) {
        double q2_thr = thr * thr;
        if (q2_thr > q2_min && q2_thr < q2_max) q2_set.insert(q2_thr);
    }
    return std::vector<double>(q2_set.begin(), q2_set.end());
}

// ================================================================
int main() {
    std::cout << "=== Grid vs APFEL++ comparison test ===" << std::endl;

    auto grid_def = pineapfel::load_grid_def("runcards/grid_dis.yaml");
    auto theory   = pineapfel::load_theory_card("runcards/theory.yaml");
    auto op_card  = pineapfel::load_operator_card("runcards/operator.yaml");

    // Derive channels (same logic as build_grid) so the test reference
    // can use them
    {
        double q2_max = 0;
        for (const auto &bin : grid_def.bins)
            q2_max = std::max(q2_max, bin.upper[0]);
        int nf_max = apfel::NF(std::sqrt(q2_max), theory.quark_thresholds);
        grid_def.channels = pineapfel::derive_channels(grid_def.process,
            grid_def.observable,
            grid_def.current,
            grid_def.cc_sign,
            nf_max);
    }

    pineappl_grid *grid  = pineapfel::build_grid(grid_def, theory, op_card);
    std::size_t    nbins = pineappl_grid_bin_count(grid);
    std::size_t    nords = pineappl_grid_order_count(grid);
    std::size_t    nchan = grid_def.channels.size();

    // Build independent APFEL++ objects
    std::vector<apfel::SubGrid> sgs;
    for (const auto &sg : op_card.xgrid)
        sgs.emplace_back(sg.n_knots, sg.x_min, sg.poly_degree);
    const apfel::Grid g{sgs};

    auto sf_init = apfel::InitializeF2NCObjectsZM(g, theory.quark_thresholds);

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

    auto   q2_nodes  = derive_q2_nodes(grid_def.bins, theory.quark_thresholds);
    bool   timelike  = (grid_def.process == pineapfel::ProcessType::SIA);

    int    failures  = 0;
    double tolerance = 1e-6;

    // ============================================================
    // TEST 1: PineAPPL LO convolution vs APFEL++ operator convolution
    //
    // Convolve the PineAPPL grid (LO only) with a toy PDF and compare
    // against independently computing the same observable:
    //   - For each Q² node: build coefficient function operators,
    //     convolve with the toy PDF distributions, sum over channels
    //   - Sum the per-Q²-node results (matching PineAPPL's sum)
    // ============================================================
    std::printf("\n--- TEST 1: PineAPPL LO convolution vs APFEL++ ---\n");
    {
        auto lo_mask = std::make_unique<bool[]>(nords);
        for (std::size_t i = 0; i < nords; i++) lo_mask[i] = (i == 0);

        std::vector<double> pineappl_lo(nbins, 0.0);
        pineappl_grid_convolve_with_one(grid,
            2212,
            xfx_callback,
            alphas_callback,
            static_cast<void *>(&as_tab),
            lo_mask.get(),
            nullptr,
            1.0,
            1.0,
            pineappl_lo.data());

        for (std::size_t ibin = 0; ibin < nbins; ibin++) {
            double x_center = std::sqrt(grid_def.bins[ibin].lower.back() *
                                        grid_def.bins[ibin].upper.back());
            double ref      = 0;

            for (double q2 : q2_nodes) {
                double Q       = std::sqrt(q2);
                int    nf      = apfel::NF(Q, theory.quark_thresholds);
                auto   charges = apfel::ElectroWeakCharges(Q, timelike);
                auto   FObjQ   = sf_init(Q, charges);
                if (FObjQ.C0.count(1) == 0) continue;
                auto ops = FObjQ.C0.at(1).GetObjects();

                for (std::size_t ich = 0; ich < nchan; ich++) {
                    apfel::Operator C_ch =
                        build_channel_operator(grid_def.channels[ich],
                            ops,
                            charges,
                            nf);

                    const auto         &ch = grid_def.channels[ich];
                    apfel::Distribution pdf_ch(g,
                        [&](double const &z) -> double {
                            double sum = 0;
                            for (std::size_t ic = 0;
                                 ic < ch.pid_combinations.size();
                                 ic++) {
                                double f_val = 1.0;
                                for (int pid : ch.pid_combinations[ic])
                                    f_val *= toy_f(pid, z);
                                sum += ch.factors[ic] * f_val;
                            }
                            return sum;
                        });

                    ref += (C_ch * pdf_ch).Evaluate(x_center);
                }
            }

            double rel_diff =
                std::abs(ref) > 1e-30
                    ? std::abs(pineappl_lo[ibin] - ref) / std::abs(ref)
                    : 0.0;

            bool ok = (rel_diff < tolerance);
            std::printf(
                "  bin %zu: pineappl=%.6e  apfel++=%.6e  rel_diff=%.2e %s\n",
                ibin,
                pineappl_lo[ibin],
                ref,
                rel_diff,
                ok ? "OK" : "FAIL");
            if (!ok) failures++;
        }
    }

    // ============================================================
    // TEST 2: PineAPPL NLO convolution vs APFEL++ operator convolution
    //
    // Same as TEST 1 but including both LO and NLO contributions.
    // PineAPPL multiplies NLO subgrids.
    // ============================================================
    std::printf("\n--- TEST 2: PineAPPL LO+NLO convolution vs APFEL++ ---\n");
    {
        auto nlo_mask = std::make_unique<bool[]>(nords);
        for (std::size_t i = 0; i < nords; i++) nlo_mask[i] = (i <= 1);

        std::vector<double> pineappl_nlo(nbins, 0.0);
        pineappl_grid_convolve_with_one(grid,
            2212,
            xfx_callback,
            alphas_callback,
            static_cast<void *>(&as_tab),
            nlo_mask.get(),
            nullptr,
            1.0,
            1.0,
            pineappl_nlo.data());

        for (std::size_t ibin = 0; ibin < nbins; ibin++) {
            double x_center = std::sqrt(grid_def.bins[ibin].lower.back() *
                                        grid_def.bins[ibin].upper.back());
            double ref      = 0;

            for (double q2 : q2_nodes) {
                double Q       = std::sqrt(q2);
                int    nf      = apfel::NF(Q, theory.quark_thresholds);
                auto   charges = apfel::ElectroWeakCharges(Q, timelike);
                auto   FObjQ   = sf_init(Q, charges);
                double as_val  = as_tab.Evaluate(Q);

                // LO contribution
                if (FObjQ.C0.count(1)) {
                    auto ops0 = FObjQ.C0.at(1).GetObjects();
                    for (std::size_t ich = 0; ich < nchan; ich++) {
                        apfel::Operator C_ch =
                            build_channel_operator(grid_def.channels[ich],
                                ops0,
                                charges,
                                nf);
                        const auto         &ch = grid_def.channels[ich];
                        apfel::Distribution pdf_ch(g,
                            [&](double const &z) -> double {
                                double sum = 0;
                                for (std::size_t ic = 0;
                                     ic < ch.pid_combinations.size();
                                     ic++) {
                                    double f_val = 1.0;
                                    for (int pid : ch.pid_combinations[ic])
                                        f_val *= toy_f(pid, z);
                                    sum += ch.factors[ic] * f_val;
                                }
                                return sum;
                            });
                        ref += (C_ch * pdf_ch).Evaluate(x_center);
                    }
                }

                // NLO contribution (multiplied by alpha_s)
                if (FObjQ.C1.count(1)) {
                    auto ops1 = FObjQ.C1.at(1).GetObjects();
                    for (std::size_t ich = 0; ich < nchan; ich++) {
                        apfel::Operator C_ch =
                            build_channel_operator(grid_def.channels[ich],
                                ops1,
                                charges,
                                nf);
                        const auto         &ch = grid_def.channels[ich];
                        apfel::Distribution pdf_ch(g,
                            [&](double const &z) -> double {
                                double sum = 0;
                                for (std::size_t ic = 0;
                                     ic < ch.pid_combinations.size();
                                     ic++) {
                                    double f_val = 1.0;
                                    for (int pid : ch.pid_combinations[ic])
                                        f_val *= toy_f(pid, z);
                                    sum += ch.factors[ic] * f_val;
                                }
                                return sum;
                            });
                        ref += as_val * (C_ch * pdf_ch).Evaluate(x_center);
                    }
                }
            }

            double rel_diff =
                std::abs(ref) > 1e-30
                    ? std::abs(pineappl_nlo[ibin] - ref) / std::abs(ref)
                    : 0.0;

            bool ok = (rel_diff < tolerance);
            std::printf(
                "  bin %zu: pineappl=%.6e  apfel++=%.6e  rel_diff=%.2e %s\n",
                ibin,
                pineappl_nlo[ibin],
                ref,
                rel_diff,
                ok ? "OK" : "FAIL");
            if (!ok) failures++;
        }
    }

    // ============================================================
    // TEST 3: Full NNLO convolution (all orders)
    // ============================================================
    std::printf(
        "\n--- TEST 3: PineAPPL full NNLO convolution vs APFEL++ ---\n");
    {
        std::vector<double> pineappl_full(nbins, 0.0);
        pineappl_grid_convolve_with_one(grid,
            2212,
            xfx_callback,
            alphas_callback,
            static_cast<void *>(&as_tab),
            nullptr,
            nullptr,
            1.0,
            1.0,
            pineappl_full.data());

        for (std::size_t ibin = 0; ibin < nbins; ibin++) {
            double x_center = std::sqrt(grid_def.bins[ibin].lower.back() *
                                        grid_def.bins[ibin].upper.back());
            double ref      = 0;

            for (double q2 : q2_nodes) {
                double Q       = std::sqrt(q2);
                int    nf      = apfel::NF(Q, theory.quark_thresholds);
                auto   charges = apfel::ElectroWeakCharges(Q, timelike);
                auto   FObjQ   = sf_init(Q, charges);
                double as_val  = as_tab.Evaluate(Q);

                // Order 0 (LO), 1 (NLO), 2 (NNLO) with alpha_s^n factors
                struct {
                    int                                         order;
                    double                                      factor;
                    std::map<int, apfel::Set<apfel::Operator>> *Cn;
                } order_info[3];

                // We need to handle the maps carefully
                for (int iord = 0; iord < 3; iord++) {
                    double as_power = std::pow(as_val, iord);
                    std::map<int, apfel::Operator> ops_map;

                    if (iord == 0 && FObjQ.C0.count(1))
                        ops_map = FObjQ.C0.at(1).GetObjects();
                    else if (iord == 1 && FObjQ.C1.count(1))
                        ops_map = FObjQ.C1.at(1).GetObjects();
                    else if (iord == 2 && FObjQ.C2.count(1))
                        ops_map = FObjQ.C2.at(1).GetObjects();
                    else continue;

                    for (std::size_t ich = 0; ich < nchan; ich++) {
                        apfel::Operator C_ch =
                            build_channel_operator(grid_def.channels[ich],
                                ops_map,
                                charges,
                                nf);
                        const auto         &ch = grid_def.channels[ich];
                        apfel::Distribution pdf_ch(g,
                            [&](double const &z) -> double {
                                double sum = 0;
                                for (std::size_t ic = 0;
                                     ic < ch.pid_combinations.size();
                                     ic++) {
                                    double f_val = 1.0;
                                    for (int pid : ch.pid_combinations[ic])
                                        f_val *= toy_f(pid, z);
                                    sum += ch.factors[ic] * f_val;
                                }
                                return sum;
                            });
                        ref += as_power * (C_ch * pdf_ch).Evaluate(x_center);
                    }
                }
            }

            double rel_diff =
                std::abs(ref) > 1e-30
                    ? std::abs(pineappl_full[ibin] - ref) / std::abs(ref)
                    : 0.0;

            bool ok = (rel_diff < tolerance);
            std::printf(
                "  bin %zu: pineappl=%.6e  apfel++=%.6e  rel_diff=%.2e %s\n",
                ibin,
                pineappl_full[ibin],
                ref,
                rel_diff,
                ok ? "OK" : "FAIL");
            if (!ok) failures++;
        }
    }

    // ============================================================
    // TEST 4: PineAPPL full convolution vs APFEL++ BuildStructureFunctions
    //
    // Use APFEL++'s BuildStructureFunctions to compute the total F2
    // via the evolution basis, and compare against PineAPPL grid
    // convolution with auto-derived channels (all flavors + gluon).
    //
    // Note: PineAPPL multiplies coefficient functions by alpha_s^n,
    // while BuildStructureFunctions uses alpha_s/(4*pi) internally.
    // We pass 4*pi*alpha_s as the coupling to compensate.
    // ============================================================
    std::printf("\n--- TEST 4: PineAPPL vs BuildStructureFunctions ---\n");
    {
        // Build the structure function Observable via APFEL++
        // InDistFunc: evolution-basis distributions (i=0..12)
        // For our toy PDF (all quarks equal, no top):
        //   GLUON(0)  = g(x)
        //   SIGMA(1)  = 2*5*f(x) = 10*f(x)
        //   T35(11)   = 10*f(x)  (5 equal flavors minus 0 top)
        //   all others = 0
        auto dist_func =
            [](int const &i, double const &x, double const & /*Q*/) -> double {
            double f     = std::pow(x, 0.5) * std::pow(1.0 - x, 3.0);
            double g_val = std::pow(x, -0.1) * std::pow(1.0 - x, 5.0);
            if (i == 0) return g_val;     // gluon
            if (i == 1) return 10.0 * f;  // Sigma
            if (i == 11) return 10.0 * f; // T35
            return 0.0;
        };

        // Alphas function: pass 4*pi*alpha_s to match PineAPPL convention
        // (BuildStructureFunctions divides by 4*pi internally)
        auto alphas_func = [&](double const &Q) -> double {
            return apfel::FourPi * as_tab.Evaluate(Q);
        };

        auto couplings_func = [&](double const &Q) -> std::vector<double> {
            return apfel::ElectroWeakCharges(Q, timelike);
        };

        auto                F2_map   = apfel::BuildStructureFunctions(sf_init,
            dist_func,
            theory.pert_ord,
            alphas_func,
            couplings_func);

        // F2_map[0] is the total structure function Observable
        auto               &F2_total = F2_map.at(0);

        // PineAPPL full convolution (all orders)
        std::vector<double> pineappl_full(nbins, 0.0);
        pineappl_grid_convolve_with_one(grid,
            2212,
            xfx_callback,
            alphas_callback,
            static_cast<void *>(&as_tab),
            nullptr,
            nullptr,
            1.0,
            1.0,
            pineappl_full.data());

        for (std::size_t ibin = 0; ibin < nbins; ibin++) {
            double x_center = std::sqrt(grid_def.bins[ibin].lower.back() *
                                        grid_def.bins[ibin].upper.back());
            double ref      = 0;

            for (double q2 : q2_nodes) {
                double Q  = std::sqrt(q2);
                ref      += F2_total.Evaluate(x_center, Q);
            }

            double rel_diff =
                std::abs(ref) > 1e-30
                    ? std::abs(pineappl_full[ibin] - ref) / std::abs(ref)
                    : 0.0;

            // Looser tolerance: PineAPPL uses interpolated subgrids while
            // BuildStructureFunctions computes exact convolutions
            double bsf_tol = 1e-2;
            bool   ok      = (rel_diff < bsf_tol);
            std::printf(
                "  bin %zu: pineappl=%.6e  BSF=%.6e  rel_diff=%.2e %s\n",
                ibin,
                pineappl_full[ibin],
                ref,
                rel_diff,
                ok ? "OK" : "FAIL");
            if (!ok) failures++;
        }
    }

    // ============================================================
    // TEST 5: CC F2 Plus — PineAPPL vs BuildStructureFunctions
    //
    // Build a CC grid and compare against APFEL++ BuildStructureFunctions
    // using the CC Plus initializer and CKM weights.
    // ============================================================
    std::printf(
        "\n--- TEST 5: CC F2 Plus PineAPPL vs BuildStructureFunctions ---\n");
    {
        auto cc_grid_def =
            pineapfel::load_grid_def("runcards/grid_dis_cc.yaml");
        auto cc_theory = pineapfel::load_theory_card("runcards/theory.yaml");

        pineappl_grid *cc_grid =
            pineapfel::build_grid(cc_grid_def, cc_theory, op_card);
        std::size_t cc_nbins = pineappl_grid_bin_count(cc_grid);

        // Build CC structure function via APFEL++ BuildStructureFunctions
        auto        cc_sf_init =
            apfel::InitializeF2CCPlusObjectsZM(g, cc_theory.quark_thresholds);

        auto cc_dist_func =
            [](int const &i, double const &x, double const & /*Q*/) -> double {
            double f     = std::pow(x, 0.5) * std::pow(1.0 - x, 3.0);
            double g_val = std::pow(x, -0.1) * std::pow(1.0 - x, 5.0);
            if (i == 0) return g_val;     // gluon
            if (i == 1) return 10.0 * f;  // Sigma
            if (i == 11) return 10.0 * f; // T35
            return 0.0;
        };

        auto cc_alphas_func = [&](double const &Q) -> double {
            return apfel::FourPi * as_tab.Evaluate(Q);
        };

        auto cc_couplings_func = [&](double const &Q) -> std::vector<double> {
            (void)Q;
            return cc_theory.ckm;
        };

        auto  F2CC_map   = apfel::BuildStructureFunctions(cc_sf_init,
            cc_dist_func,
            cc_theory.pert_ord,
            cc_alphas_func,
            cc_couplings_func);

        auto &F2CC_total = F2CC_map.at(0);

        // Derive Q2 nodes and channels for the CC grid
        {
            double q2max = 0;
            for (const auto &bin : cc_grid_def.bins)
                q2max = std::max(q2max, bin.upper[0]);
            int nfm = apfel::NF(std::sqrt(q2max), cc_theory.quark_thresholds);
            cc_grid_def.channels =
                pineapfel::derive_channels(cc_grid_def.process,
                    cc_grid_def.observable,
                    cc_grid_def.current,
                    cc_grid_def.cc_sign,
                    nfm);
        }

        auto cc_q2_nodes =
            derive_q2_nodes(cc_grid_def.bins, cc_theory.quark_thresholds);

        std::vector<double> pineappl_cc(cc_nbins, 0.0);
        pineappl_grid_convolve_with_one(cc_grid,
            2212,
            xfx_callback,
            alphas_callback,
            static_cast<void *>(&as_tab),
            nullptr,
            nullptr,
            1.0,
            1.0,
            pineappl_cc.data());

        for (std::size_t ibin = 0; ibin < cc_nbins; ibin++) {
            double x_center = std::sqrt(cc_grid_def.bins[ibin].lower.back() *
                                        cc_grid_def.bins[ibin].upper.back());
            double ref      = 0;

            for (double q2 : cc_q2_nodes) {
                double Q  = std::sqrt(q2);
                ref      += F2CC_total.Evaluate(x_center, Q);
            }

            double rel_diff =
                std::abs(ref) > 1e-30
                    ? std::abs(pineappl_cc[ibin] - ref) / std::abs(ref)
                    : 0.0;

            double bsf_tol = 1e-2;
            bool   ok      = (rel_diff < bsf_tol);
            std::printf(
                "  bin %zu: pineappl=%.6e  BSF=%.6e  rel_diff=%.2e %s\n",
                ibin,
                pineappl_cc[ibin],
                ref,
                rel_diff,
                ok ? "OK" : "FAIL");
            if (!ok) failures++;
        }

        pineappl_grid_delete(cc_grid);
    }

    // ============================================================
    // TEST 6: SIDIS F2 — PineAPPL LO+NLO vs APFEL++ DoubleObject
    //
    // Build a SIDIS grid and compare PineAPPL convolution (with 2
    // convolution functions: PDF ⊗ FF) against independently
    // evaluating the DoubleObject coefficient functions from
    // InitializeSIDIS (computed in a separate TU to avoid ODR
    // violations from APFEL++'s header-defined functions).
    // ============================================================
    std::printf("\n--- TEST 6: SIDIS F2 LO+NLO PineAPPL vs APFEL++ ---\n");
    {
        auto sidis_grid_def =
            pineapfel::load_grid_def("runcards/grid_sidis.yaml");
        auto sidis_theory = pineapfel::load_theory_card("runcards/theory.yaml");

        pineappl_grid *sidis_grid =
            pineapfel::build_grid(sidis_grid_def, sidis_theory, op_card);
        std::size_t sidis_nbins = pineappl_grid_bin_count(sidis_grid);

        auto        sidis_q2_nodes =
            derive_q2_nodes(sidis_grid_def.bins, sidis_theory.quark_thresholds);

        // PineAPPL convolution with 2 convolution functions (PDF and FF).
        // Both use the same toy_f functional form.
        std::vector<double> pineappl_sidis(sidis_nbins, 0.0);
        void               *pdfs_state[2] = {nullptr, nullptr};

        pineappl_grid_convolve(sidis_grid,
            xfx_callback,
            alphas_callback,
            pdfs_state,
            static_cast<void *>(&as_tab),
            nullptr, // all orders
            nullptr, // all channels
            nullptr, // no bin indices
            0,       // nb_scales
            nullptr, // mu_scales
            pineappl_sidis.data());

        // Reference computation via helper (separate TU with SIDIS.h)
        std::vector<std::vector<double>> bin_x_bounds, bin_z_bounds;
        for (const auto &bin : sidis_grid_def.bins) {
            bin_x_bounds.push_back({bin.lower[1], bin.upper[1]});
            bin_z_bounds.push_back({bin.lower[2], bin.upper[2]});
        }

        auto alphas_func = [&](double Q) -> double {
            return as_tab.Evaluate(Q);
        };

        auto ref_vals = compute_sidis_reference(g,
            sidis_theory.quark_thresholds,
            sidis_q2_nodes,
            bin_x_bounds,
            bin_z_bounds,
            {},
            1, // max_alpha_s = 1 (LO + NLO)
            toy_f,
            alphas_func);

        for (std::size_t ibin = 0; ibin < sidis_nbins; ibin++) {
            double ref = ref_vals[ibin];
            double rel_diff =
                std::abs(ref) > 1e-30
                    ? std::abs(pineappl_sidis[ibin] - ref) / std::abs(ref)
                    : 0.0;

            bool ok = (rel_diff < tolerance);
            std::printf(
                "  bin %zu: pineappl=%.6e  apfel++=%.6e  rel_diff=%.2e %s\n",
                ibin,
                pineappl_sidis[ibin],
                ref,
                rel_diff,
                ok ? "OK" : "FAIL");
            if (!ok) failures++;
        }

        pineappl_grid_delete(sidis_grid);
    }

    // ============================================================
    // TEST 7: Polarized DIS g₁ — PineAPPL vs BuildStructureFunctions
    //
    // Build a polarized DIS grid (g₁) and compare against APFEL++
    // BuildStructureFunctions using Initializeg1NCObjectsZM.
    // Follows the same pattern as TEST 4.
    // ============================================================
    std::printf(
        "\n--- TEST 7: Polarized DIS g1 PineAPPL vs BuildStructureFunctions "
        "---\n");
    {
        auto pol_grid_def =
            pineapfel::load_grid_def("runcards/grid_dis_pol.yaml");
        auto pol_theory = pineapfel::load_theory_card("runcards/theory.yaml");

        pineappl_grid *pol_grid =
            pineapfel::build_grid(pol_grid_def, pol_theory, op_card);
        std::size_t pol_nbins = pineappl_grid_bin_count(pol_grid);

        // Build g1 structure function via APFEL++ BuildStructureFunctions
        auto        g1_sf_init =
            apfel::Initializeg1NCObjectsZM(g, pol_theory.quark_thresholds);

        // Same dist_func as TEST 4 (equal toy quarks mapped to evolution basis)
        auto pol_dist_func =
            [](int const &i, double const &x, double const & /*Q*/) -> double {
            double f     = std::pow(x, 0.5) * std::pow(1.0 - x, 3.0);
            double g_val = std::pow(x, -0.1) * std::pow(1.0 - x, 5.0);
            if (i == 0) return g_val;     // gluon
            if (i == 1) return 10.0 * f;  // Sigma
            if (i == 11) return 10.0 * f; // T35
            return 0.0;
        };

        auto pol_alphas_func = [&](double const &Q) -> double {
            return apfel::FourPi * as_tab.Evaluate(Q);
        };

        auto pol_couplings_func = [&](double const &Q) -> std::vector<double> {
            return apfel::ElectroWeakCharges(Q, false);
        };

        auto                g1_map = apfel::BuildStructureFunctions(g1_sf_init,
            pol_dist_func,
            pol_theory.pert_ord,
            pol_alphas_func,
            pol_couplings_func);
        auto               &g1_total = g1_map.at(0);

        // PineAPPL full convolution (all orders).
        // Must use pineappl_grid_convolve (not convolve_with_one) because
        // the grid has PolPDF convolution type; convolve_with_one always
        // assumes UnpolPDF and would panic at runtime.
        std::vector<double> pineappl_pol(pol_nbins, 0.0);
        void               *pol_state[1] = {nullptr};

        pineappl_grid_convolve(pol_grid,
            xfx_callback,
            alphas_callback,
            pol_state,
            static_cast<void *>(&as_tab),
            nullptr, // all orders
            nullptr, // all channels
            nullptr, // no bin indices
            0,       // nb_scales
            nullptr, // mu_scales
            pineappl_pol.data());

        auto pol_q2_nodes =
            derive_q2_nodes(pol_grid_def.bins, pol_theory.quark_thresholds);

        for (std::size_t ibin = 0; ibin < pol_nbins; ibin++) {
            double x_center = std::sqrt(pol_grid_def.bins[ibin].lower.back() *
                                        pol_grid_def.bins[ibin].upper.back());
            double ref      = 0;

            for (double q2 : pol_q2_nodes) {
                double Q  = std::sqrt(q2);
                ref      += g1_total.Evaluate(x_center, Q);
            }

            double rel_diff =
                std::abs(ref) > 1e-30
                    ? std::abs(pineappl_pol[ibin] - ref) / std::abs(ref)
                    : 0.0;

            double bsf_tol = 1e-2;
            bool   ok      = (rel_diff < bsf_tol);
            std::printf(
                "  bin %zu: pineappl=%.6e  BSF=%.6e  rel_diff=%.2e %s\n",
                ibin,
                pineappl_pol[ibin],
                ref,
                rel_diff,
                ok ? "OK" : "FAIL");
            if (!ok) failures++;
        }

        pineappl_grid_delete(pol_grid);
    }

    // ============================================================
    // TEST 8: Polarized SIDIS G₁ — PineAPPL LO+NLO vs APFEL++ DoubleObject
    //
    // Build a polarized SIDIS grid (G₁) and compare PineAPPL
    // convolution (POL_PDF ⊗ UNPOL_FF) against independently
    // evaluating the DoubleObject coefficient functions from
    // InitializeSIDISpol.  Follows the same pattern as TEST 6.
    // ============================================================
    std::printf(
        "\n--- TEST 8: Polarized SIDIS G1 LO+NLO PineAPPL vs APFEL++ ---\n");
    {
        auto pol_sidis_grid_def =
            pineapfel::load_grid_def("runcards/grid_sidis_pol.yaml");
        auto pol_sidis_theory =
            pineapfel::load_theory_card("runcards/theory.yaml");

        pineappl_grid *pol_sidis_grid =
            pineapfel::build_grid(pol_sidis_grid_def,
                pol_sidis_theory,
                op_card);
        std::size_t pol_sidis_nbins = pineappl_grid_bin_count(pol_sidis_grid);

        auto pol_sidis_q2_nodes     = derive_q2_nodes(pol_sidis_grid_def.bins,
            pol_sidis_theory.quark_thresholds);

        // PineAPPL convolution with 2 convolution functions (POL_PDF ⊗
        // UNPOL_FF)
        std::vector<double> pineappl_pol_sidis(pol_sidis_nbins, 0.0);
        void               *pdfs_state[2] = {nullptr, nullptr};

        pineappl_grid_convolve(pol_sidis_grid,
            xfx_callback,
            alphas_callback,
            pdfs_state,
            static_cast<void *>(&as_tab),
            nullptr, // all orders
            nullptr, // all channels
            nullptr, // no bin indices
            0,       // nb_scales
            nullptr, // mu_scales
            pineappl_pol_sidis.data());

        // Reference computation via helper (separate TU with SIDISpol.h)
        std::vector<std::vector<double>> pol_bin_x_bounds, pol_bin_z_bounds;
        for (const auto &bin : pol_sidis_grid_def.bins) {
            pol_bin_x_bounds.push_back({bin.lower[1], bin.upper[1]});
            pol_bin_z_bounds.push_back({bin.lower[2], bin.upper[2]});
        }

        auto pol_alphas_func = [&](double Q) -> double {
            return as_tab.Evaluate(Q);
        };

        auto pol_ref_vals = compute_sidis_pol_reference(g,
            pol_sidis_theory.quark_thresholds,
            pol_sidis_q2_nodes,
            pol_bin_x_bounds,
            pol_bin_z_bounds,
            1, // max_alpha_s = 1 (LO + NLO)
            toy_f,
            pol_alphas_func);

        for (std::size_t ibin = 0; ibin < pol_sidis_nbins; ibin++) {
            double ref = pol_ref_vals[ibin];
            double rel_diff =
                std::abs(ref) > 1e-30
                    ? std::abs(pineappl_pol_sidis[ibin] - ref) / std::abs(ref)
                    : 0.0;

            bool ok = (rel_diff < tolerance);
            std::printf(
                "  bin %zu: pineappl=%.6e  apfel++=%.6e  rel_diff=%.2e %s\n",
                ibin,
                pineappl_pol_sidis[ibin],
                ref,
                rel_diff,
                ok ? "OK" : "FAIL");
            if (!ok) failures++;
        }

        pineappl_grid_delete(pol_sidis_grid);
    }

    // ============================================================
    // TEST 9: FFN F2 — PineAPPL vs BuildStructureFunctions (Massive)
    //
    // Build a DIS grid with MassScheme: FFN and compare against
    // APFEL++ BuildStructureFunctions using InitializeF2NCObjectsMassive.
    // ============================================================
    std::printf(
        "\n--- TEST 9: FFN F2 PineAPPL vs BuildStructureFunctions ---\n");
    {
        auto ffn_grid_def =
            pineapfel::load_grid_def("runcards/grid_dis_ffn.yaml");
        auto ffn_theory = pineapfel::load_theory_card("runcards/theory.yaml");

        pineappl_grid *ffn_grid =
            pineapfel::build_grid(ffn_grid_def, ffn_theory, op_card);
        std::size_t ffn_nbins   = pineappl_grid_bin_count(ffn_grid);

        auto        ffn_sf_init = apfel::InitializeF2NCObjectsMassive(g,
            ffn_theory.heavy_quark_masses);

        auto        ffn_dist_func =
            [](int const &i, double const &x, double const & /*Q*/) -> double {
            double f     = std::pow(x, 0.5) * std::pow(1.0 - x, 3.0);
            double g_val = std::pow(x, -0.1) * std::pow(1.0 - x, 5.0);
            if (i == 0) return g_val;
            if (i == 1) return 10.0 * f;
            if (i == 11) return 10.0 * f;
            return 0.0;
        };

        auto ffn_alphas_func = [&](double const &Q) -> double {
            return apfel::FourPi * as_tab.Evaluate(Q);
        };

        auto ffn_couplings_func = [&](double const &Q) -> std::vector<double> {
            return apfel::ElectroWeakCharges(Q, false);
        };

        auto  F2FFN_map   = apfel::BuildStructureFunctions(ffn_sf_init,
            ffn_dist_func,
            ffn_theory.pert_ord,
            ffn_alphas_func,
            ffn_couplings_func);
        auto &F2FFN_total = F2FFN_map.at(0);

        std::vector<double> pineappl_ffn(ffn_nbins, 0.0);
        void               *ffn_state[1] = {nullptr};
        pineappl_grid_convolve(ffn_grid,
            xfx_callback,
            alphas_callback,
            ffn_state,
            static_cast<void *>(&as_tab),
            nullptr,
            nullptr,
            nullptr,
            0,
            nullptr,
            pineappl_ffn.data());

        auto ffn_q2_nodes =
            derive_q2_nodes(ffn_grid_def.bins, ffn_theory.quark_thresholds);

        for (std::size_t ibin = 0; ibin < ffn_nbins; ibin++) {
            double x_center = std::sqrt(ffn_grid_def.bins[ibin].lower.back() *
                                        ffn_grid_def.bins[ibin].upper.back());
            double ref      = 0;
            for (double q2 : ffn_q2_nodes) {
                double Q  = std::sqrt(q2);
                ref      += F2FFN_total.Evaluate(x_center, Q);
            }
            double rel_diff =
                std::abs(ref) > 1e-30
                    ? std::abs(pineappl_ffn[ibin] - ref) / std::abs(ref)
                    : 0.0;
            double bsf_tol = 1e-2;
            bool   ok      = (rel_diff < bsf_tol);
            std::printf(
                "  bin %zu: pineappl=%.6e  BSF=%.6e  rel_diff=%.2e %s\n",
                ibin,
                pineappl_ffn[ibin],
                ref,
                rel_diff,
                ok ? "OK" : "FAIL");
            if (!ok) failures++;
        }
        pineappl_grid_delete(ffn_grid);
    }

    // ============================================================
    // TEST 10: MassiveZero F2 — PineAPPL vs BuildStructureFunctions
    //
    // Build a DIS grid with MassScheme: MassiveZero and compare against
    // APFEL++ BuildStructureFunctions using InitializeF2NCObjectsMassiveZero.
    // ============================================================
    std::printf("\n--- TEST 10: MassiveZero F2 PineAPPL vs "
                "BuildStructureFunctions ---\n");
    {
        auto mz_grid_def =
            pineapfel::load_grid_def("runcards/grid_dis_mz.yaml");
        auto mz_theory = pineapfel::load_theory_card("runcards/theory.yaml");

        pineappl_grid *mz_grid =
            pineapfel::build_grid(mz_grid_def, mz_theory, op_card);
        std::size_t mz_nbins   = pineappl_grid_bin_count(mz_grid);

        auto        mz_sf_init = apfel::InitializeF2NCObjectsMassiveZero(g,
            mz_theory.heavy_quark_masses);

        auto        mz_dist_func =
            [](int const &i, double const &x, double const & /*Q*/) -> double {
            double f     = std::pow(x, 0.5) * std::pow(1.0 - x, 3.0);
            double g_val = std::pow(x, -0.1) * std::pow(1.0 - x, 5.0);
            if (i == 0) return g_val;
            if (i == 1) return 10.0 * f;
            if (i == 11) return 10.0 * f;
            return 0.0;
        };

        auto mz_alphas_func = [&](double const &Q) -> double {
            return apfel::FourPi * as_tab.Evaluate(Q);
        };

        auto mz_couplings_func = [&](double const &Q) -> std::vector<double> {
            return apfel::ElectroWeakCharges(Q, false);
        };

        auto  F2MZ_map   = apfel::BuildStructureFunctions(mz_sf_init,
            mz_dist_func,
            mz_theory.pert_ord,
            mz_alphas_func,
            mz_couplings_func);
        auto &F2MZ_total = F2MZ_map.at(0);

        std::vector<double> pineappl_mz(mz_nbins, 0.0);
        void               *mz_state[1] = {nullptr};
        pineappl_grid_convolve(mz_grid,
            xfx_callback,
            alphas_callback,
            mz_state,
            static_cast<void *>(&as_tab),
            nullptr,
            nullptr,
            nullptr,
            0,
            nullptr,
            pineappl_mz.data());

        auto mz_q2_nodes =
            derive_q2_nodes(mz_grid_def.bins, mz_theory.quark_thresholds);

        for (std::size_t ibin = 0; ibin < mz_nbins; ibin++) {
            double x_center = std::sqrt(mz_grid_def.bins[ibin].lower.back() *
                                        mz_grid_def.bins[ibin].upper.back());
            double ref      = 0;
            for (double q2 : mz_q2_nodes) {
                double Q  = std::sqrt(q2);
                ref      += F2MZ_total.Evaluate(x_center, Q);
            }
            double rel_diff =
                std::abs(ref) > 1e-30
                    ? std::abs(pineappl_mz[ibin] - ref) / std::abs(ref)
                    : 0.0;
            double bsf_tol = 1e-2;
            bool   ok      = (rel_diff < bsf_tol);
            std::printf(
                "  bin %zu: pineappl=%.6e  BSF=%.6e  rel_diff=%.2e %s\n",
                ibin,
                pineappl_mz[ibin],
                ref,
                rel_diff,
                ok ? "OK" : "FAIL");
            if (!ok) failures++;
        }
        pineappl_grid_delete(mz_grid);
    }

    // ============================================================
    // TEST 11: FONLL F2 — PineAPPL vs manual combination
    //
    // Build a DIS grid with MassScheme: FONLL and compare against
    // the manually combined reference: F_FONLL = F_ZM + F_FFN - F_MZ.
    // ============================================================
    std::printf("\n--- TEST 11: FONLL F2 PineAPPL vs manual ZM+FFN-MZ ---\n");
    {
        auto fonll_grid_def =
            pineapfel::load_grid_def("runcards/grid_dis_fonll.yaml");
        auto fonll_theory = pineapfel::load_theory_card("runcards/theory.yaml");

        pineappl_grid *fonll_grid =
            pineapfel::build_grid(fonll_grid_def, fonll_theory, op_card);
        std::size_t fonll_nbins    = pineappl_grid_bin_count(fonll_grid);

        // ZM initializer (already built as sf_init above)
        // FFN initializer
        auto        fonll_ffn_init = apfel::InitializeF2NCObjectsMassive(g,
            fonll_theory.heavy_quark_masses);
        // MassiveZero initializer
        auto        fonll_mz_init  = apfel::InitializeF2NCObjectsMassiveZero(g,
            fonll_theory.heavy_quark_masses);

        auto        fonll_dist_func =
            [](int const &i, double const &x, double const & /*Q*/) -> double {
            double f     = std::pow(x, 0.5) * std::pow(1.0 - x, 3.0);
            double g_val = std::pow(x, -0.1) * std::pow(1.0 - x, 5.0);
            if (i == 0) return g_val;
            if (i == 1) return 10.0 * f;
            if (i == 11) return 10.0 * f;
            return 0.0;
        };

        auto fonll_alphas_func = [&](double const &Q) -> double {
            return apfel::FourPi * as_tab.Evaluate(Q);
        };

        auto fonll_couplings_func =
            [&](double const &Q) -> std::vector<double> {
            return apfel::ElectroWeakCharges(Q, false);
        };

        auto  F2ZM_map_f  = apfel::BuildStructureFunctions(sf_init,
            fonll_dist_func,
            fonll_theory.pert_ord,
            fonll_alphas_func,
            fonll_couplings_func);
        auto  F2FFN_map_f = apfel::BuildStructureFunctions(fonll_ffn_init,
            fonll_dist_func,
            fonll_theory.pert_ord,
            fonll_alphas_func,
            fonll_couplings_func);
        auto  F2MZ_map_f  = apfel::BuildStructureFunctions(fonll_mz_init,
            fonll_dist_func,
            fonll_theory.pert_ord,
            fonll_alphas_func,
            fonll_couplings_func);

        auto &F2ZM_tot    = F2ZM_map_f.at(0);
        auto &F2FFN_tot   = F2FFN_map_f.at(0);
        auto &F2MZ_tot    = F2MZ_map_f.at(0);

        std::vector<double> pineappl_fonll(fonll_nbins, 0.0);
        void               *fonll_state[1] = {nullptr};
        pineappl_grid_convolve(fonll_grid,
            xfx_callback,
            alphas_callback,
            fonll_state,
            static_cast<void *>(&as_tab),
            nullptr,
            nullptr,
            nullptr,
            0,
            nullptr,
            pineappl_fonll.data());

        auto fonll_q2_nodes =
            derive_q2_nodes(fonll_grid_def.bins, fonll_theory.quark_thresholds);

        for (std::size_t ibin = 0; ibin < fonll_nbins; ibin++) {
            double x_center = std::sqrt(fonll_grid_def.bins[ibin].lower.back() *
                                        fonll_grid_def.bins[ibin].upper.back());
            double ref      = 0;
            for (double q2 : fonll_q2_nodes) {
                double Q  = std::sqrt(q2);
                // F_FONLL = F_ZM + F_FFN - F_MZ
                ref      += F2ZM_tot.Evaluate(x_center, Q) +
                       F2FFN_tot.Evaluate(x_center, Q) -
                       F2MZ_tot.Evaluate(x_center, Q);
            }
            double rel_diff =
                std::abs(ref) > 1e-30
                    ? std::abs(pineappl_fonll[ibin] - ref) / std::abs(ref)
                    : 0.0;
            double bsf_tol = 1e-2;
            bool   ok      = (rel_diff < bsf_tol);
            std::printf(
                "  bin %zu: pineappl=%.6e  ref=%.6e  rel_diff=%.2e %s\n",
                ibin,
                pineappl_fonll[ibin],
                ref,
                rel_diff,
                ok ? "OK" : "FAIL");
            if (!ok) failures++;
        }
        pineappl_grid_delete(fonll_grid);
    }

    // ============================================================
    // Summary
    // ============================================================
    std::printf("\n=== Summary: %d failures ===\n", failures);
    pineappl_grid_delete(grid);

    if (failures > 0) {
        std::printf("TEST FAILED\n");
        return 1;
    }
    std::printf("TEST PASSED\n");
    return 0;
}
