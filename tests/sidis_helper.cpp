#include "sidis_helper.h"
#include <pineapfel/grid_gen.h>
#include <pineapfel/sidis_api.h>

#include <cmath>
#include <utility>

// Same 8-point Gauss–Legendre rule as pineapfel/src/fill.cpp for
// ∫_{z_lo}^{z_hi} of the external SIDIS z variable.
static constexpr int    kSidisZQuadN                = 8;
static constexpr double kSidisZQuadXi[kSidisZQuadN] = {
    -0.960289856497536231683560853,
    -0.796666477413626739591553936,
    -0.525532409916328985817739049,
    -0.183434642495649804939476142,
    0.183434642495649804939476142,
    0.525532409916328985817739049,
    0.796666477413626739591553936,
    0.960289856497536231683560853,
};
static constexpr double kSidisZQuadWi[kSidisZQuadN] = {
    0.1012285362903762591506828,
    0.2223810343749745229420239,
    0.3137066459079072956537718,
    0.3626837833783619786074925,
    0.3626837833783619786074925,
    0.3137066459079072956537718,
    0.2223810343749745229420239,
    0.1012285362903762591506828,
};

// APFEL++ sidisbuilder uses fixed electromagnetic charge factors QCh2
// for SIDIS channels (not ElectroWeakCharges(Q,false)).
static constexpr double kSidisQCh2[6] =
    {4. / 9., 1. / 9., 1. / 9., 4. / 9., 1. / 9., 4. / 9.};

static void sidis_z_quadrature_point(double z_lo,
    double                                  z_hi,
    int                                     k,
    double                                 &z_out,
    double                                 &w_out) {
    const double mid  = 0.5 * (z_hi + z_lo);
    const double half = 0.5 * (z_hi - z_lo);
    z_out             = mid + half * kSidisZQuadXi[k];
    w_out             = half * kSidisZQuadWi[k];
}

// Helper: evaluate a DoubleOperator contracted with a PDF and an FF at fixed
// (x, z) using direct kernel access (avoids DistributionOperator::Evaluate
// which crashes on null matrix data).
static double eval_dop_at(const apfel::Grid &g,
    const apfel::LagrangeInterpolator       &li,
    const apfel::DoubleOperator             &dop,
    int                                      pid_pdf,
    int                                      pid_ff,
    double                                   x,
    double                                   z,
    std::function<double(int, double)>       toy_f) {
    const apfel::Grid &g1  = dop.GetFirstGrid();
    const apfel::Grid &g2  = dop.GetSecondGrid();
    const auto        &K   = dop.GetDoubleOperator();

    const int          ng1 = g1.nGrids();
    const int          ng2 = g2.nGrids();

    int                nx  = 0;
    for (int ig = 0; ig < ng1; ig++) nx += g1.GetSubGrid(ig).nx();
    int nz = 0;
    for (int ig = 0; ig < ng2; ig++) nz += g2.GetSubGrid(ig).nx();

    const auto &x_nodes = g.GetJointGrid().GetGrid();

    int         ig1_xc  = ng1 - 1;
    for (int ig = 0; ig < ng1; ig++)
        if (g1.GetSubGrid(ig).xMin() > x) {
            ig1_xc = ig - 1;
            break;
        }

    int ig2_zc = ng2 - 1;
    for (int ig = 0; ig < ng2; ig++)
        if (g2.GetSubGrid(ig).xMin() > z) {
            ig2_zc = ig - 1;
            break;
        }

    const apfel::SubGrid &sg1    = g1.GetSubGrid(ig1_xc);
    const apfel::SubGrid &sg2    = g2.GetSubGrid(ig2_zc);
    const int             nx_sub = sg1.nx();
    const int             nz_sub = sg2.nx();

    const auto            bx     = li.SumBounds(x, sg1);
    const int             lo_x = bx[0], hi_x = bx[1];
    std::vector<double>   W_x(hi_x - lo_x);
    for (int b = 0; b < hi_x - lo_x; b++)
        W_x[b] = li.Interpolant(lo_x + b, x, sg1);

    const auto          bz   = li.SumBounds(z, sg2);
    const int           lo_z = bz[0], hi_z = bz[1];
    std::vector<double> W_z(hi_z - lo_z);
    for (int b = 0; b < hi_z - lo_z; b++)
        W_z[b] = li.Interpolant(lo_z + b, z, sg2);

    int x_jg_off = 0;
    for (int ig = 0; ig < ig1_xc; ig++) x_jg_off += g1.GetSubGrid(ig).nx();
    int z_jg_off = 0;
    for (int ig = 0; ig < ig2_zc; ig++) z_jg_off += g2.GetSubGrid(ig).nx();

    const auto &K_sub  = K[ig1_xc][ig2_zc];
    double      result = 0.0;

    for (int b_ox = 0; b_ox < hi_x - lo_x; b_ox++) {
        const int    out_x = lo_x + b_ox;
        const double wx    = W_x[b_ox];
        for (int b_oz = 0; b_oz < hi_z - lo_z; b_oz++) {
            const int    out_z = lo_z + b_oz;
            const double wxz   = wx * W_z[b_oz];
            for (int alpha = 0; alpha <= nx_sub - out_x - 1; alpha++) {
                const int inp_x = x_jg_off + out_x + alpha;
                if (inp_x >= nx) break;
                // Skip APFL++ Lagrange extension nodes beyond x = 1.
                // fill.cpp drops these from the PineAPPL subgrid; the
                // reference must apply the same cut to remain consistent.
                if (x_nodes[inp_x] > 1.0) continue;
                for (int gamma = 0; gamma <= nz_sub - out_z - 1; gamma++) {
                    const int inp_z = z_jg_off + out_z + gamma;
                    if (inp_z >= nz) break;
                    if (x_nodes[inp_z] > 1.0) continue;
                    result += wxz * K_sub(0, alpha)(0, gamma) * x_nodes[inp_x] *
                              x_nodes[inp_z] * toy_f(pid_pdf, x_nodes[inp_x]) *
                              toy_f(pid_ff, x_nodes[inp_z]);
                }
            }
        }
    }

    return result;
}

std::vector<double> compute_sidis_reference(const apfel::Grid &g,
    const apfel::SidisNNLOObjects                             &sobj,
    const std::vector<double>                                 &thresholds,
    const std::vector<double>                                 &q2_nodes,
    const std::vector<std::vector<double>>                    &bin_x_bounds,
    const std::vector<std::vector<double>>                    &bin_z_bounds,
    const std::vector<double> & /*charges_dummy*/,
    int                                max_alpha_s,
    std::function<double(int, double)> toy_f,
    std::function<double(double)>      alphas_func) {
    size_t                            nbins = bin_x_bounds.size();

    const apfel::LagrangeInterpolator li{g};

    double                            q2_max_val = 0;
    for (double q2 : q2_nodes) q2_max_val = std::max(q2_max_val, q2);
    int  nf_max   = apfel::NF(std::sqrt(q2_max_val), thresholds);

    auto channels = pineapfel::derive_channels(pineapfel::ProcessType::SIDIS,
        pineapfel::Observable::F2,
        pineapfel::Current::NC,
        pineapfel::CCSign::Plus,
        nf_max);

    // Key lookup helpers for FT operators (LO/NLO only here)
    auto ft_key   = [](int alpha_s, int channel_type) -> std::string {
        if (alpha_s == 0 && channel_type == 0) return "DoubleIdentity";
        if (alpha_s == 1 && channel_type == 0) return "C1TQ2Q";
        if (alpha_s == 1 && channel_type == 1) return "C1TQ2G";
        if (alpha_s == 1 && channel_type == 2) return "C1TG2Q";
        return "";
    };

    std::vector<double> result(nbins, 0.0);

    for (size_t ibin = 0; ibin < nbins; ibin++) {
        double x_center =
            std::sqrt(bin_x_bounds[ibin][0] * bin_x_bounds[ibin][1]);
        const double z_lo = bin_z_bounds[ibin][0];
        const double z_hi = bin_z_bounds[ibin][1];

        for (double q2 : q2_nodes) {
            double Q      = std::sqrt(q2);
            int    nf     = apfel::NF(Q, thresholds);
            double as_val = alphas_func(Q);

            for (int izq = 0; izq < kSidisZQuadN; izq++) {
                double z_ext = 0.0, dz_w = 0.0;
                sidis_z_quadrature_point(z_lo, z_hi, izq, z_ext, dz_w);
                const auto &x_nodes = g.GetJointGrid().GetGrid();
                if (z_ext < x_nodes.front() || z_ext > 1.0) continue;

                for (size_t ich = 0; ich < channels.size(); ich++) {
                    int quark_idx    = static_cast<int>(ich / 3);
                    int channel_type = static_cast<int>(ich % 3);

                    if (quark_idx + 1 > nf) continue;
                    double      e_q_sq = kSidisQCh2[quark_idx];
                    const auto &ch     = channels[ich];

                    for (int alpha_s = 0; alpha_s <= max_alpha_s; alpha_s++) {
                        std::string key = ft_key(alpha_s, channel_type);
                        if (key.empty()) continue;
                        auto it = sobj.FT.find(key);
                        if (it == sobj.FT.end()) continue;

                        double as_power =
                            std::pow(as_val / (4.0 * M_PI), alpha_s);
                        for (size_t ic = 0; ic < ch.pid_combinations.size();
                             ic++) {
                            int pid_pdf   = ch.pid_combinations[ic][0];
                            int pid_ff    = ch.pid_combinations[ic][1];
                            result[ibin] += dz_w * as_power * e_q_sq *
                                            ch.factors[ic] /
                                            (x_center * z_ext) *
                                            eval_dop_at(g,
                                                li,
                                                it->second,
                                                pid_pdf,
                                                pid_ff,
                                                x_center,
                                                z_ext,
                                                toy_f);
                        }
                    }
                }
            }
        }
    }

    return result;
}

std::vector<double> compute_sidis_pol_reference(const apfel::Grid &g,
    const apfel::SidisNNLOObjects                                 &sobj,
    const std::vector<double>                                     &thresholds,
    const std::vector<double>                                     &q2_nodes,
    const std::vector<std::vector<double>>                        &bin_x_bounds,
    const std::vector<std::vector<double>>                        &bin_z_bounds,
    int                                                            max_alpha_s,
    std::function<double(int, double)>                             toy_f,
    std::function<double(double)> alphas_func) {
    size_t                            nbins = bin_x_bounds.size();

    const apfel::LagrangeInterpolator li{g};

    double                            q2_max_val = 0;
    for (double q2 : q2_nodes) q2_max_val = std::max(q2_max_val, q2);
    int  nf_max   = apfel::NF(std::sqrt(q2_max_val), thresholds);

    auto channels = pineapfel::derive_channels(pineapfel::ProcessType::SIDIS,
        pineapfel::Observable::F2,
        pineapfel::Current::NC,
        pineapfel::CCSign::Plus,
        nf_max);

    // Key lookup helpers for G1 operators (LO/NLO only here)
    auto g1_key   = [](int alpha_s, int channel_type) -> std::string {
        if (alpha_s == 0 && channel_type == 0) return "DoubleIdentity";
        if (alpha_s == 1 && channel_type == 0) return "DC1Q2Q";
        if (alpha_s == 1 && channel_type == 1) return "DC1Q2G";
        if (alpha_s == 1 && channel_type == 2) return "DC1G2Q";
        return "";
    };

    std::vector<double> result(nbins, 0.0);

    for (size_t ibin = 0; ibin < nbins; ibin++) {
        double x_center =
            std::sqrt(bin_x_bounds[ibin][0] * bin_x_bounds[ibin][1]);
        const double z_lo = bin_z_bounds[ibin][0];
        const double z_hi = bin_z_bounds[ibin][1];

        for (double q2 : q2_nodes) {
            double Q      = std::sqrt(q2);
            int    nf     = apfel::NF(Q, thresholds);
            double as_val = alphas_func(Q);

            for (int izq = 0; izq < kSidisZQuadN; izq++) {
                double z_ext = 0.0, dz_w = 0.0;
                sidis_z_quadrature_point(z_lo, z_hi, izq, z_ext, dz_w);
                const auto &x_nodes = g.GetJointGrid().GetGrid();
                if (z_ext < x_nodes.front() || z_ext > 1.0) continue;

                for (size_t ich = 0; ich < channels.size(); ich++) {
                    int quark_idx    = static_cast<int>(ich / 3);
                    int channel_type = static_cast<int>(ich % 3);

                    if (quark_idx + 1 > nf) continue;
                    double      e_q_sq = kSidisQCh2[quark_idx];
                    const auto &ch     = channels[ich];

                    for (int alpha_s = 0; alpha_s <= max_alpha_s; alpha_s++) {
                        std::string key = g1_key(alpha_s, channel_type);
                        if (key.empty()) continue;
                        auto it = sobj.G1.find(key);
                        if (it == sobj.G1.end()) continue;

                        double as_power =
                            std::pow(as_val / (4.0 * M_PI), alpha_s);
                        for (size_t ic = 0; ic < ch.pid_combinations.size();
                             ic++) {
                            int pid_pdf   = ch.pid_combinations[ic][0];
                            int pid_ff    = ch.pid_combinations[ic][1];
                            result[ibin] += dz_w * as_power * e_q_sq *
                                            ch.factors[ic] /
                                            (x_center * z_ext) *
                                            eval_dop_at(g,
                                                li,
                                                it->second,
                                                pid_pdf,
                                                pid_ff,
                                                x_center,
                                                z_ext,
                                                toy_f);
                        }
                    }
                }
            }
        }
    }

    return result;
}

// ============================================================
// NNLO SIDIS F2 reference via direct DoubleOperator evaluation.
// ============================================================
std::vector<double> compute_sidis_nnlo_reference(const apfel::Grid &g,
    const apfel::SidisNNLOObjects                                  &sobj,
    const std::vector<double>                                      &thresholds,
    int                                                             nf_max,
    const std::vector<double>                                      &q2_nodes,
    const std::vector<std::vector<double>> &bin_x_bounds,
    const std::vector<std::vector<double>> &bin_z_bounds,
    std::function<double(int, double)>      toy_f,
    std::function<double(double)>           alphas_func) {
    constexpr double QCh[6] =
        {2. / 3., -1. / 3., -1. / 3., 2. / 3., -1. / 3., 2. / 3.};

    size_t                            nbins = bin_x_bounds.size();
    std::vector<double>               result(nbins, 0.0);

    const apfel::LagrangeInterpolator li{g};

    for (size_t ibin = 0; ibin < nbins; ibin++) {
        double x_c = std::sqrt(bin_x_bounds[ibin][0] * bin_x_bounds[ibin][1]);
        const double z_lo = bin_z_bounds[ibin][0];
        const double z_hi = bin_z_bounds[ibin][1];

        for (double q2 : q2_nodes) {
            double Q      = std::sqrt(q2);
            int    nf     = apfel::NF(Q, thresholds);
            double as_val = alphas_func(Q);

            double sumch2 = 0;
            for (int q = 0; q < nf; q++) sumch2 += kSidisQCh2[q];

            const std::string snf = "_nf" + std::to_string(nf);

            for (int izq = 0; izq < kSidisZQuadN; izq++) {
                double z_ext = 0.0, dz_w = 0.0;
                sidis_z_quadrature_point(z_lo, z_hi, izq, z_ext, dz_w);
                const auto &jg = g.GetJointGrid().GetGrid();
                if (z_ext < jg.front() || z_ext > 1.0) continue;

                double as2 = dz_w * (as_val / (4.0 * M_PI)) *
                             (as_val / (4.0 * M_PI)) / (x_c * z_ext);

                auto eval_dop = [&](const apfel::DoubleOperator &dop,
                                    int                          pid_pdf,
                                    int pid_ff) -> double {
                    return eval_dop_at(g,
                        li,
                        dop,
                        pid_pdf,
                        pid_ff,
                        x_c,
                        z_ext,
                        toy_f);
                };

                using Combo     = std::pair<int, int>;
                auto sum_combos = [&](const apfel::DoubleOperator &dop,
                                      const std::vector<Combo>    &combos,
                                      const std::vector<double>   &factors) {
                    double s = 0;
                    for (size_t ic = 0; ic < combos.size(); ic++)
                        s += factors[ic] *
                             eval_dop(dop, combos[ic].first, combos[ic].second);
                    return s;
                };

                auto find_ft = [&](const std::string &key)
                    -> const apfel::DoubleOperator * {
                    auto it = sobj.FT.find(key);
                    return it != sobj.FT.end() ? &it->second : nullptr;
                };

                // ns/qq per quark
                if (const auto *op = find_ft("C2TQ2QNS" + snf))
                    for (int q = 1; q <= nf_max; q++) {
                        if (q > nf) continue;
                        result[ibin] += as2 * kSidisQCh2[q - 1] *
                                        sum_combos(*op,
                                            {
                                                { q,  q},
                                                {-q, -q}
                        },
                                            {1., 1.});
                    }

                // gq per quark
                if (const auto *op = find_ft("C2TQ2G"))
                    for (int q = 1; q <= nf_max; q++) {
                        if (q > nf) continue;
                        result[ibin] += as2 * kSidisQCh2[q - 1] *
                                        sum_combos(*op,
                                            {
                                                { q, 21},
                                                {-q, 21}
                        },
                                            {1., 1.});
                    }

                // qg per quark
                if (const auto *op = find_ft("C2TG2Q"))
                    for (int q = 1; q <= nf_max; q++) {
                        if (q > nf) continue;
                        result[ibin] += as2 * kSidisQCh2[q - 1] *
                                        sum_combos(*op,
                                            {
                                                {21,  q},
                                                {21, -q}
                        },
                                            {1., 1.});
                    }

                // gg
                if (const auto *op = find_ft("C2TG2G"))
                    result[ibin] += as2 * sumch2 *
                                    sum_combos(*op,
                                        {
                                            {21, 21}
                    },
                                        {1.});

                // ps per quark
                if (const auto *op = find_ft("C2TQ2QPS"))
                    for (int q = 1; q <= nf_max; q++) {
                        if (q > nf) continue;
                        result[ibin] += as2 * sumch2 *
                                        sum_combos(*op,
                                            {
                                                { q,  q},
                                                {-q, -q}
                        },
                                            {1., 1.});
                    }

                // qbq per quark
                if (const auto *op = find_ft("C2TQ2QB"))
                    for (int q = 1; q <= nf_max; q++) {
                        if (q > nf) continue;
                        result[ibin] += as2 * kSidisQCh2[q - 1] *
                                        sum_combos(*op,
                                            {
                                                { q, -q},
                                                {-q,  q}
                        },
                                            {1., 1.});
                    }

                // qpq1 per source j
                if (const auto *op = find_ft("C2TQ2QP1"))
                    if (nf_max >= 2)
                        for (int j = 1; j <= nf_max; j++) {
                            if (j > nf) continue;
                            std::vector<Combo>  combos;
                            std::vector<double> factors;
                            for (int k = 1; k <= nf_max; k++) {
                                if (k == j) continue;
                                combos.insert(combos.end(),
                                    {
                                        { j,  k},
                                        { j, -k},
                                        {-j,  k},
                                        {-j, -k}
                                });
                                factors.insert(factors.end(), {1., 1., 1., 1.});
                            }
                            result[ibin] += as2 * kSidisQCh2[j - 1] *
                                            sum_combos(*op, combos, factors);
                        }

                // qpq2 per target i
                if (const auto *op = find_ft("C2TQ2QP2"))
                    if (nf_max >= 2)
                        for (int i = 1; i <= nf_max; i++) {
                            if (i > nf) continue;
                            std::vector<Combo>  combos;
                            std::vector<double> factors;
                            for (int k = 1; k <= nf_max; k++) {
                                if (k == i) continue;
                                combos.insert(combos.end(),
                                    {
                                        { k,  i},
                                        { k, -i},
                                        {-k,  i},
                                        {-k, -i}
                                });
                                factors.insert(factors.end(), {1., 1., 1., 1.});
                            }
                            result[ibin] += as2 * kSidisQCh2[i - 1] *
                                            sum_combos(*op, combos, factors);
                        }

                // qpq3 per ordered pair (a,b) with a!=b
                if (const auto *op = find_ft("C2TQ2QP3"))
                    for (int a = 1; a <= nf_max; a++) {
                        if (a > nf) continue;
                        for (int b = 1; b <= nf_max; b++) {
                            if (b > nf || b == a) continue;
                            result[ibin] += as2 * QCh[a - 1] * QCh[b - 1] *
                                            sum_combos(*op,
                                                {
                                                    { a,  b},
                                                    {-a,  b},
                                                    { a, -b},
                                                    {-a, -b}
                            },
                                                {1., -1., -1., 1.});
                        }
                    }
            }
        }
    }

    return result;
}

// ============================================================
// NNLO SIDIS G1 (polarised) — same layout as compute_sidis_nnlo_reference.
// ============================================================
std::vector<double> compute_sidis_g1_nnlo_reference(const apfel::Grid &g,
    const apfel::SidisNNLOObjects                                     &sobj,
    const std::vector<double>              &thresholds,
    int                                     nf_max,
    const std::vector<double>              &q2_nodes,
    const std::vector<std::vector<double>> &bin_x_bounds,
    const std::vector<std::vector<double>> &bin_z_bounds,
    std::function<double(int, double)>      toy_f,
    std::function<double(double)>           alphas_func) {
    constexpr double QCh[6] =
        {2. / 3., -1. / 3., -1. / 3., 2. / 3., -1. / 3., 2. / 3.};

    size_t                            nbins = bin_x_bounds.size();
    std::vector<double>               result(nbins, 0.0);

    const apfel::LagrangeInterpolator li{g};

    for (size_t ibin = 0; ibin < nbins; ibin++) {
        double x_c = std::sqrt(bin_x_bounds[ibin][0] * bin_x_bounds[ibin][1]);
        const double z_lo = bin_z_bounds[ibin][0];
        const double z_hi = bin_z_bounds[ibin][1];

        for (double q2 : q2_nodes) {
            double Q      = std::sqrt(q2);
            int    nf     = apfel::NF(Q, thresholds);
            double as_val = alphas_func(Q);

            double sumch2 = 0;
            for (int q = 0; q < nf; q++) sumch2 += kSidisQCh2[q];

            const std::string snf = "_nf" + std::to_string(nf);

            for (int izq = 0; izq < kSidisZQuadN; izq++) {
                double z_ext = 0.0, dz_w = 0.0;
                sidis_z_quadrature_point(z_lo, z_hi, izq, z_ext, dz_w);
                const auto &jg = g.GetJointGrid().GetGrid();
                if (z_ext < jg.front() || z_ext > 1.0) continue;

                double as2 = dz_w * (as_val / (4.0 * M_PI)) *
                             (as_val / (4.0 * M_PI)) / (x_c * z_ext);

                auto eval_dop = [&](const apfel::DoubleOperator &dop,
                                    int                          pid_pdf,
                                    int pid_ff) -> double {
                    return eval_dop_at(g,
                        li,
                        dop,
                        pid_pdf,
                        pid_ff,
                        x_c,
                        z_ext,
                        toy_f);
                };

                using Combo     = std::pair<int, int>;
                auto sum_combos = [&](const apfel::DoubleOperator &dop,
                                      const std::vector<Combo>    &combos,
                                      const std::vector<double>   &factors) {
                    double s = 0;
                    for (size_t ic = 0; ic < combos.size(); ic++)
                        s += factors[ic] *
                             eval_dop(dop, combos[ic].first, combos[ic].second);
                    return s;
                };

                auto find_g1 = [&](const std::string &key)
                    -> const apfel::DoubleOperator * {
                    auto it = sobj.G1.find(key);
                    return it != sobj.G1.end() ? &it->second : nullptr;
                };

                if (const auto *op = find_g1("DC2Q2QNS" + snf))
                    for (int q = 1; q <= nf_max; q++) {
                        if (q > nf) continue;
                        result[ibin] += as2 * kSidisQCh2[q - 1] *
                                        sum_combos(*op,
                                            {
                                                { q,  q},
                                                {-q, -q}
                        },
                                            {1., 1.});
                    }

                if (const auto *op = find_g1("DC2Q2G" + snf))
                    for (int q = 1; q <= nf_max; q++) {
                        if (q > nf) continue;
                        result[ibin] += as2 * kSidisQCh2[q - 1] *
                                        sum_combos(*op,
                                            {
                                                { q, 21},
                                                {-q, 21}
                        },
                                            {1., 1.});
                    }

                if (const auto *op = find_g1("DC2G2Q" + snf))
                    for (int q = 1; q <= nf_max; q++) {
                        if (q > nf) continue;
                        result[ibin] += as2 * kSidisQCh2[q - 1] *
                                        sum_combos(*op,
                                            {
                                                {21,  q},
                                                {21, -q}
                        },
                                            {1., 1.});
                    }

                if (const auto *op = find_g1("DC2G2G"))
                    result[ibin] += as2 * sumch2 *
                                    sum_combos(*op,
                                        {
                                            {21, 21}
                    },
                                        {1.});

                if (const auto *op = find_g1("DC2Q2QPS"))
                    for (int q = 1; q <= nf_max; q++) {
                        if (q > nf) continue;
                        result[ibin] += as2 * sumch2 *
                                        sum_combos(*op,
                                            {
                                                { q,  q},
                                                {-q, -q}
                        },
                                            {1., 1.});
                    }

                if (const auto *op = find_g1("DC2Q2QB"))
                    for (int q = 1; q <= nf_max; q++) {
                        if (q > nf) continue;
                        result[ibin] += as2 * kSidisQCh2[q - 1] *
                                        sum_combos(*op,
                                            {
                                                { q, -q},
                                                {-q,  q}
                        },
                                            {1., 1.});
                    }

                if (const auto *op = find_g1("DC2Q2QP1"))
                    if (nf_max >= 2)
                        for (int j = 1; j <= nf_max; j++) {
                            if (j > nf) continue;
                            std::vector<Combo>  combos;
                            std::vector<double> factors;
                            for (int k = 1; k <= nf_max; k++) {
                                if (k == j) continue;
                                combos.insert(combos.end(),
                                    {
                                        { j,  k},
                                        { j, -k},
                                        {-j,  k},
                                        {-j, -k}
                                });
                                factors.insert(factors.end(), {1., 1., 1., 1.});
                            }
                            result[ibin] += as2 * kSidisQCh2[j - 1] *
                                            sum_combos(*op, combos, factors);
                        }

                if (const auto *op = find_g1("DC2Q2QP2"))
                    if (nf_max >= 2)
                        for (int i = 1; i <= nf_max; i++) {
                            if (i > nf) continue;
                            std::vector<Combo>  combos;
                            std::vector<double> factors;
                            for (int k = 1; k <= nf_max; k++) {
                                if (k == i) continue;
                                combos.insert(combos.end(),
                                    {
                                        { k,  i},
                                        { k, -i},
                                        {-k,  i},
                                        {-k, -i}
                                });
                                factors.insert(factors.end(), {1., 1., 1., 1.});
                            }
                            result[ibin] += as2 * kSidisQCh2[i - 1] *
                                            sum_combos(*op, combos, factors);
                        }

                if (const auto *op = find_g1("DC2Q2QP3"))
                    for (int a = 1; a <= nf_max; a++) {
                        if (a > nf) continue;
                        for (int b = 1; b <= nf_max; b++) {
                            if (b > nf || b == a) continue;
                            result[ibin] += as2 * QCh[a - 1] * QCh[b - 1] *
                                            sum_combos(*op,
                                                {
                                                    { a,  b},
                                                    {-a,  b},
                                                    { a, -b},
                                                    {-a, -b}
                            },
                                                {1., -1., -1., 1.});
                        }
                    }
            }
        }
    }

    return result;
}
