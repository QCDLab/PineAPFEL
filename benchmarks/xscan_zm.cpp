// Benchmark: PineAPFEL vs APFEL++ x-scan at fixed Q²=10 (ZM-VFN).
//
// Builds PineAPFEL grids for DIS, SIA, and SIDIS (NC, F2, ZM-VFN scheme) with
// 50 narrow log-spaced bins in x (or z for SIA) at Q²=[10.0, 10.001], then
// compares the convolved PineAPPL result against the APFEL++ reference obtained
// via direct operator convolution (same method as fill.cpp / sidis_helper).
//
// Output: benchmarks/xscan_zm.dat

#include "sidis_helper.h"
#include <apfel/apfelxx.h>
#include <pineapfel.h>
#include <pineappl_capi.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>

// TODO: Replace with a proper LHAPDF/NeoPDF call.
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
    return static_cast<apfel::TabulateObject<double> *>(state)->Evaluate(
        std::sqrt(q2));
}

// Build combined channel operator from APFEL++ NC DIS operator objects.
// Mirrors `build_channel_operator` in test_grid_vs_apfelxx.cpp.
static apfel::Operator build_channel_operator(
    const pineapfel::ChannelDef          &channel,
    const std::map<int, apfel::Operator> &ops,
    const std::vector<double>            &charges,
    int                                   nf) {
    const apfel::Operator &CNS    = ops.at(apfel::DISNCBasis::CNS);
    const apfel::Operator &CS     = ops.at(apfel::DISNCBasis::CS);
    const apfel::Operator &CG     = ops.at(apfel::DISNCBasis::CG);

    double                 sum_ch = 0;
    for (int i = 0; i < nf && i < (int)charges.size(); i++)
        sum_ch += charges[i];

    bool is_gluon = false;
    for (const auto &combo : channel.pid_combinations)
        if (combo.size() == 1 && combo[0] == 21) {
            is_gluon = true;
            break;
        }

    if (is_gluon) return sum_ch * CG;

    int quark_pid = 0;
    for (const auto &combo : channel.pid_combinations)
        if (combo.size() == 1 && combo[0] > 0) {
            quark_pid = combo[0];
            break;
        }
    if (quark_pid == 0)
        for (const auto &combo : channel.pid_combinations)
            if (combo.size() == 1 && combo[0] < 0) {
                quark_pid = -combo[0];
                break;
            }

    int    q_idx = quark_pid - 1;
    double e_q_sq =
        (q_idx >= 0 && q_idx < (int)charges.size()) ? charges[q_idx] : 0.0;
    return e_q_sq * CNS + (sum_ch / 6.0) * (CS - CNS);
}

// Reproduce `derive_q2_nodes` from fill.cpp
static std::vector<double> derive_q2_nodes(
    const std::vector<pineapfel::BinDef> &bins,
    const std::vector<double>            &thresholds,
    int                                   n_intermediate = 3) {
    std::set<double> q2_set;
    for (const auto &bin : bins) {
        double lo = bin.lower[0], hi = bin.upper[0];
        q2_set.insert(lo);
        q2_set.insert(hi);
        if (n_intermediate > 0) {
            double lll = std::log(lo), llh = std::log(hi);
            for (int i = 1; i <= n_intermediate; i++) {
                double f = static_cast<double>(i) / (n_intermediate + 1);
                q2_set.insert(std::exp(lll + f * (llh - lll)));
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

static double rel_diff(double pineappl, double ref) {
    return std::abs(ref) > 1e-30 ? std::abs(pineappl - ref) / std::abs(ref)
                                 : 0.0;
}

int main() {
    auto theory  = pineapfel::load_theory_card("runcards/theory.yaml");
    auto op_card = pineapfel::load_operator_card("runcards/operator.yaml");

    std::vector<apfel::SubGrid> sgs;
    for (const auto &sg : op_card.xgrid)
        sgs.emplace_back(sg.n_knots, sg.x_min, sg.poly_degree);
    const apfel::Grid             g{sgs};

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

    constexpr double              Q2_lo = 10.0, Q2_hi = 10.001;

    constexpr int                 NX = 50;
    constexpr double X_MIN = 1e-6, X_MAX = 1.0; // DIS / SIDIS x range
    constexpr double Z_MIN = 1e-2, Z_MAX = 1.0; // SIA z range

    auto             make_edges = [](int n, double lo, double hi) {
        std::vector<double> e(n + 1);
        for (int i = 0; i <= n; i++)
            e[i] = std::exp(std::log(lo) + static_cast<double>(i) / n *
                                               (std::log(hi) - std::log(lo)));
        return e;
    };
    const auto x_edges      = make_edges(NX, X_MIN, X_MAX);
    const auto z_edges      = make_edges(NX, Z_MIN, Z_MAX);

    // Build 1-D bins [Q2_lo,Q2_hi]×[e_i, e_{i+1}]
    auto       make_bins_1d = [&](const std::vector<double> &edges) {
        std::vector<pineapfel::BinDef> bins;
        for (int i = 0; i < NX; i++)
            bins.push_back({
                {Q2_lo,     edges[i]},
                {Q2_hi, edges[i + 1]}
            });
        return bins;
    };

    // SIDIS: 3-D bins [Q2_lo,Q2_hi]×[x_i,x_{i+1}]×[z_lo,z_hi]
    constexpr double SIDIS_Z_LO = 0.005, SIDIS_Z_HI = 0.015;
    auto             make_bins_sidis = [&]() {
        std::vector<pineapfel::BinDef> bins;
        for (int i = 0; i < NX; i++)
            bins.push_back({
                {Q2_lo,     x_edges[i], SIDIS_Z_LO},
                {Q2_hi, x_edges[i + 1], SIDIS_Z_HI}
            });
        return bins;
    };

    const std::vector<double>              norms(NX, 1.0);
    const std::vector<pineapfel::OrderDef> orders_nnlo = {
        {0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0},
        {2, 0, 0, 0, 0}
    };
    const std::vector<pineapfel::OrderDef> orders_nlo = {
        {0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0}
    };

    std::ofstream out("benchmarks/xscan_zm.dat");
    out << "# process  x_or_z  pineappl  apfelxx  rel_diff\n";
    out << std::scientific;

    // =========================================================================
    // DIS NC F2 ZM — x-scan
    // =========================================================================
    std::cout << "\n=== DIS NC F2 ZM — x-scan ===" << std::endl;
    {
        pineapfel::GridDef gd;
        gd.process           = pineapfel::ProcessType::DIS;
        gd.observable        = pineapfel::Observable::F2;
        gd.current           = pineapfel::Current::NC;
        gd.mass_scheme       = pineapfel::MassScheme::ZM;
        gd.pid_basis         = PINEAPPL_PID_BASIS_PDG;
        gd.hadron_pids       = {2212};
        gd.convolution_types = {PINEAPPL_CONV_TYPE_UNPOL_PDF};
        gd.orders            = orders_nnlo;
        gd.bins              = make_bins_1d(x_edges);
        gd.normalizations    = norms;

        pineappl_grid *grid  = pineapfel::build_grid(gd, theory, op_card);

        // Derive channels (mirrors build_grid internal logic)
        {
            double q2_max = 0;
            for (const auto &bin : gd.bins)
                q2_max = std::max(q2_max, bin.upper[0]);
            int nf_max  = apfel::NF(std::sqrt(q2_max), theory.quark_thresholds);
            gd.channels = pineapfel::derive_channels(gd.process,
                gd.observable,
                gd.current,
                gd.cc_sign,
                nf_max);
        }
        std::size_t nchan = gd.channels.size();

        // APFEL++ reference via direct operator convolution (same as fill.cpp)
        auto        zm_sf_init =
            apfel::InitializeF2NCObjectsZM(g, theory.quark_thresholds);

        auto q2_nodes = derive_q2_nodes(gd.bins, theory.quark_thresholds);

        std::vector<double> result(NX, 0.0);
        pineappl_grid_convolve_with_one(grid,
            2212,
            xfx_callback,
            alphas_callback,
            static_cast<void *>(&as_tab),
            nullptr,
            nullptr,
            1.0,
            1.0,
            result.data());

        for (int i = 0; i < NX; i++) {
            double xc  = std::sqrt(x_edges[i] * x_edges[i + 1]);
            double ref = 0.0;
            for (double q2 : q2_nodes) {
                double Q       = std::sqrt(q2);
                int    nf      = apfel::NF(Q, theory.quark_thresholds);
                auto   charges = apfel::ElectroWeakCharges(Q, false);
                auto   FObjQ   = zm_sf_init(Q, charges);
                double as_val  = as_tab.Evaluate(Q);
                for (int iord = 0; iord < 3; iord++) {
                    const auto &C = (iord == 0)   ? FObjQ.C0
                                    : (iord == 1) ? FObjQ.C1
                                                  : FObjQ.C2;
                    if (C.count(1) == 0) continue;
                    double as_power = std::pow(as_val, iord);
                    auto   ops      = C.at(1).GetObjects();
                    for (std::size_t ich = 0; ich < nchan; ich++) {
                        apfel::Operator C_ch =
                            build_channel_operator(gd.channels[ich],
                                ops,
                                charges,
                                nf);
                        const auto         &ch = gd.channels[ich];
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
                        ref += as_power * (C_ch * pdf_ch).Evaluate(xc);
                    }
                }
            }
            double rd = rel_diff(result[i], ref);
            std::printf("  x=%.4e  pineappl=%.6e  apfel++=%.6e  rd=%.2e\n",
                xc,
                result[i],
                ref,
                rd);
            out << "DIS  " << xc << "  " << result[i] << "  " << ref << "  "
                << rd << "\n";
        }
        pineappl_grid_delete(grid);
    }

    // =========================================================================
    // SIA NC F2 ZM — z-scan
    // =========================================================================
    std::cout << "\n=== SIA NC F2 ZM — z-scan ===" << std::endl;
    {
        pineapfel::GridDef gd;
        gd.process           = pineapfel::ProcessType::SIA;
        gd.observable        = pineapfel::Observable::F2;
        gd.current           = pineapfel::Current::NC;
        gd.mass_scheme       = pineapfel::MassScheme::ZM;
        gd.pid_basis         = PINEAPPL_PID_BASIS_PDG;
        gd.hadron_pids       = {211};
        gd.convolution_types = {PINEAPPL_CONV_TYPE_UNPOL_FF};
        gd.orders            = orders_nnlo;
        gd.bins              = make_bins_1d(z_edges);
        gd.normalizations    = norms;

        pineappl_grid *grid  = pineapfel::build_grid(gd, theory, op_card);

        // Derive channels (mirrors build_grid internal logic)
        {
            double q2_max = 0;
            for (const auto &bin : gd.bins)
                q2_max = std::max(q2_max, bin.upper[0]);
            int nf_max  = apfel::NF(std::sqrt(q2_max), theory.quark_thresholds);
            gd.channels = pineapfel::derive_channels(gd.process,
                gd.observable,
                gd.current,
                gd.cc_sign,
                nf_max);
        }
        std::size_t nchan = gd.channels.size();

        // APFEL++ reference via direct operator convolution (timelike ZM)
        auto        sia_sf_init =
            apfel::InitializeF2NCObjectsZMT(g, theory.quark_thresholds);

        auto q2_nodes = derive_q2_nodes(gd.bins, theory.quark_thresholds);

        std::vector<double> result(NX, 0.0);
        void               *sia_state[1] = {nullptr};
        pineappl_grid_convolve(grid,
            xfx_callback,
            alphas_callback,
            sia_state,
            static_cast<void *>(&as_tab),
            nullptr,
            nullptr,
            nullptr,
            0,
            nullptr,
            result.data());

        for (int i = 0; i < NX; i++) {
            double zc  = std::sqrt(z_edges[i] * z_edges[i + 1]);
            double ref = 0.0;
            for (double q2 : q2_nodes) {
                double Q       = std::sqrt(q2);
                int    nf      = apfel::NF(Q, theory.quark_thresholds);
                auto   charges = apfel::ElectroWeakCharges(Q, true); // timelike
                auto   FObjQ   = sia_sf_init(Q, charges);
                double as_val  = as_tab.Evaluate(Q);
                for (int iord = 0; iord < 3; iord++) {
                    const auto &C = (iord == 0)   ? FObjQ.C0
                                    : (iord == 1) ? FObjQ.C1
                                                  : FObjQ.C2;
                    if (C.count(1) == 0) continue;
                    double as_power = std::pow(as_val, iord);
                    auto   ops      = C.at(1).GetObjects();
                    for (std::size_t ich = 0; ich < nchan; ich++) {
                        apfel::Operator C_ch =
                            build_channel_operator(gd.channels[ich],
                                ops,
                                charges,
                                nf);
                        const auto         &ch = gd.channels[ich];
                        apfel::Distribution ff_ch(g,
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
                        ref += as_power * (C_ch * ff_ch).Evaluate(zc);
                    }
                }
            }
            double rd = rel_diff(result[i], ref);
            std::printf("  z=%.4e  pineappl=%.6e  apfel++=%.6e  rd=%.2e\n",
                zc,
                result[i],
                ref,
                rd);
            out << "SIA  " << zc << "  " << result[i] << "  " << ref << "  "
                << rd << "\n";
        }
        pineappl_grid_delete(grid);
    }

    // =========================================================================
    // SIDIS NC F2 ZM — x-scan at fixed z=[0.005, 0.015]
    // =========================================================================
    std::cout << "\n=== SIDIS NC F2 ZM — x-scan ===" << std::endl;
    {
        pineapfel::GridDef gd;
        gd.process           = pineapfel::ProcessType::SIDIS;
        gd.observable        = pineapfel::Observable::F2;
        gd.current           = pineapfel::Current::NC;
        gd.mass_scheme       = pineapfel::MassScheme::ZM;
        gd.pid_basis         = PINEAPPL_PID_BASIS_PDG;
        gd.hadron_pids       = {2212, 211};
        gd.convolution_types = {PINEAPPL_CONV_TYPE_UNPOL_PDF,
            PINEAPPL_CONV_TYPE_UNPOL_FF};
        gd.orders            = orders_nlo;
        gd.bins              = make_bins_sidis();
        gd.normalizations    = norms;

        pineappl_grid *grid  = pineapfel::build_grid(gd, theory, op_card);

        // APFEL++ reference via sidis_helper
        auto q2_nodes = derive_q2_nodes(gd.bins, theory.quark_thresholds);
        std::vector<std::vector<double>> bx(NX), bz(NX);
        for (int i = 0; i < NX; i++) {
            bx[i] = {x_edges[i], x_edges[i + 1]};
            bz[i] = {SIDIS_Z_LO, SIDIS_Z_HI};
        }
        auto alphas_func = [&](double Q) { return as_tab.Evaluate(Q); };
        auto ref_vals    = compute_sidis_reference(g,
            theory.quark_thresholds,
            q2_nodes,
            bx,
            bz,
               {},
            1 /*LO+NLO*/,
            toy_f,
            alphas_func);

        std::vector<double> result(NX, 0.0);
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
            result.data());

        for (int i = 0; i < NX; i++) {
            double xc  = std::sqrt(x_edges[i] * x_edges[i + 1]);
            double ref = ref_vals[i];
            double rd  = rel_diff(result[i], ref);
            std::printf("  x=%.4e  pineappl=%.6e  apfel++=%.6e  rd=%.2e\n",
                xc,
                result[i],
                ref,
                rd);
            out << "SIDIS  " << xc << "  " << result[i] << "  " << ref << "  "
                << rd << "\n";
        }
        pineappl_grid_delete(grid);
    }

    out.close();
    std::cout << "\nResults written to benchmarks/xscan_zm.dat" << std::endl;
    return 0;
}
