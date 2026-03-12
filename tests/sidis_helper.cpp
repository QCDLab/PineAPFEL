#include "sidis_helper.h"
#include <pineapfel/grid_gen.h>
#include <pineapfel/sidis_api.h>

#include <cmath>
#include <utility>

// Helper: evaluate a DoubleOperator contracted with a PDF and an FF at fixed
// (x, z). Computes  (dop · MultiplyFirstBy(pdf)).Evaluate(x) · ff  evaluated at
// z.
static double eval_dop_at(const apfel::Grid &g,
    const apfel::DoubleOperator             &dop,
    int                                      pid_pdf,
    int                                      pid_ff,
    double                                   x,
    double                                   z,
    std::function<double(int, double)>       toy_f) {
    apfel::Distribution pdf_d(g,
        [&](double const &xx) { return toy_f(pid_pdf, xx); });
    apfel::Distribution ff_d(g,
        [&](double const &zz) { return toy_f(pid_ff, zz); });
    return (dop.MultiplyFirstBy(pdf_d).Evaluate(x) * ff_d).Evaluate(z);
}

std::vector<double> compute_sidis_reference(const apfel::Grid &g,
    const std::vector<double>                                 &thresholds,
    const std::vector<double>                                 &q2_nodes,
    const std::vector<std::vector<double>>                    &bin_x_bounds,
    const std::vector<std::vector<double>>                    &bin_z_bounds,
    const std::vector<double> & /*charges_dummy*/,
    int                                max_alpha_s,
    std::function<double(int, double)> toy_f,
    std::function<double(double)>      alphas_func) {
    auto   sobj       = pineapfel::init_sidis_nnlo(g, thresholds);
    size_t nbins      = bin_x_bounds.size();

    double q2_max_val = 0;
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
        double z_center =
            std::sqrt(bin_z_bounds[ibin][0] * bin_z_bounds[ibin][1]);

        for (double q2 : q2_nodes) {
            double Q       = std::sqrt(q2);
            int    nf      = apfel::NF(Q, thresholds);
            auto   charges = apfel::ElectroWeakCharges(Q, false);
            double as_val  = alphas_func(Q);

            for (size_t ich = 0; ich < channels.size(); ich++) {
                int quark_idx    = static_cast<int>(ich / 3);
                int channel_type = static_cast<int>(ich % 3);

                if (quark_idx + 1 > nf) continue;
                double      e_q_sq = charges[quark_idx];
                const auto &ch     = channels[ich];

                for (int alpha_s = 0; alpha_s <= max_alpha_s; alpha_s++) {
                    std::string key = ft_key(alpha_s, channel_type);
                    if (key.empty()) continue;
                    auto it = sobj.FT.find(key);
                    if (it == sobj.FT.end()) continue;

                    double as_power = std::pow(as_val, alpha_s);
                    for (size_t ic = 0; ic < ch.pid_combinations.size(); ic++) {
                        int pid_pdf   = ch.pid_combinations[ic][0];
                        int pid_ff    = ch.pid_combinations[ic][1];
                        result[ibin] += as_power * e_q_sq * ch.factors[ic] *
                                        eval_dop_at(g,
                                            it->second,
                                            pid_pdf,
                                            pid_ff,
                                            x_center,
                                            z_center,
                                            toy_f);
                    }
                }
            }
        }
    }

    return result;
}

std::vector<double> compute_sidis_pol_reference(const apfel::Grid &g,
    const std::vector<double>                                     &thresholds,
    const std::vector<double>                                     &q2_nodes,
    const std::vector<std::vector<double>>                        &bin_x_bounds,
    const std::vector<std::vector<double>>                        &bin_z_bounds,
    int                                                            max_alpha_s,
    std::function<double(int, double)>                             toy_f,
    std::function<double(double)> alphas_func) {
    auto   sobj       = pineapfel::init_sidis_nnlo(g, thresholds);
    size_t nbins      = bin_x_bounds.size();

    double q2_max_val = 0;
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
        double z_center =
            std::sqrt(bin_z_bounds[ibin][0] * bin_z_bounds[ibin][1]);

        for (double q2 : q2_nodes) {
            double Q       = std::sqrt(q2);
            int    nf      = apfel::NF(Q, thresholds);
            auto   charges = apfel::ElectroWeakCharges(Q, false);
            double as_val  = alphas_func(Q);

            for (size_t ich = 0; ich < channels.size(); ich++) {
                int quark_idx    = static_cast<int>(ich / 3);
                int channel_type = static_cast<int>(ich % 3);

                if (quark_idx + 1 > nf) continue;
                double      e_q_sq = charges[quark_idx];
                const auto &ch     = channels[ich];

                for (int alpha_s = 0; alpha_s <= max_alpha_s; alpha_s++) {
                    std::string key = g1_key(alpha_s, channel_type);
                    if (key.empty()) continue;
                    auto it = sobj.G1.find(key);
                    if (it == sobj.G1.end()) continue;

                    double as_power = std::pow(as_val, alpha_s);
                    for (size_t ic = 0; ic < ch.pid_combinations.size(); ic++) {
                        int pid_pdf   = ch.pid_combinations[ic][0];
                        int pid_ff    = ch.pid_combinations[ic][1];
                        result[ibin] += as_power * e_q_sq * ch.factors[ic] *
                                        eval_dop_at(g,
                                            it->second,
                                            pid_pdf,
                                            pid_ff,
                                            x_center,
                                            z_center,
                                            toy_f);
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

    size_t              nbins = bin_x_bounds.size();
    std::vector<double> result(nbins, 0.0);

    for (size_t ibin = 0; ibin < nbins; ibin++) {
        double x_c = std::sqrt(bin_x_bounds[ibin][0] * bin_x_bounds[ibin][1]);
        double z_c = std::sqrt(bin_z_bounds[ibin][0] * bin_z_bounds[ibin][1]);

        for (double q2 : q2_nodes) {
            double Q       = std::sqrt(q2);
            int    nf      = apfel::NF(Q, thresholds);
            auto   charges = apfel::ElectroWeakCharges(Q, false);
            double as2     = alphas_func(Q) * alphas_func(Q);

            double sumch2  = 0;
            for (int q = 0; q < nf; q++) sumch2 += charges[q];

            const std::string snf      = "_nf" + std::to_string(nf);

            auto              eval_dop = [&](const apfel::DoubleOperator &dop,
                                int                          pid_pdf,
                                int                          pid_ff) -> double {
                return eval_dop_at(g, dop, pid_pdf, pid_ff, x_c, z_c, toy_f);
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

            auto find_ft =
                [&](const std::string &key) -> const apfel::DoubleOperator * {
                auto it = sobj.FT.find(key);
                return it != sobj.FT.end() ? &it->second : nullptr;
            };

            // ns/qq per quark
            if (const auto *op = find_ft("C2TQ2QNS" + snf))
                for (int q = 1; q <= nf_max; q++) {
                    if (q > nf) continue;
                    result[ibin] += as2 * charges[q - 1] *
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
                    result[ibin] += as2 * charges[q - 1] *
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
                    result[ibin] += as2 * charges[q - 1] *
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
                    result[ibin] += as2 * charges[q - 1] *
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
                        result[ibin] += as2 * charges[j - 1] *
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
                        result[ibin] += as2 * charges[i - 1] *
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

    return result;
}
