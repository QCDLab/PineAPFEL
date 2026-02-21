#include "sidis_helper.h"
#include <grid_gen.h>
#include <sidis_api.h>

#include <cmath>

std::vector<double> compute_sidis_reference(const apfel::Grid &g,
    const std::vector<double>                                 &thresholds,
    const std::vector<double>                                 &q2_nodes,
    const std::vector<std::vector<double>>                    &bin_x_bounds,
    const std::vector<std::vector<double>>                    &bin_z_bounds,
    const std::vector<double> & /*charges_dummy*/,
    int                                max_alpha_s,
    std::function<double(int, double)> toy_f,
    std::function<double(double)>      alphas_func) {
    auto   sobj       = pineapfel::init_sidis(g, thresholds);
    size_t nbins      = bin_x_bounds.size();

    // Determine nf_max
    double q2_max_val = 0;
    for (double q2 : q2_nodes) q2_max_val = std::max(q2_max_val, q2);
    int  nf_max   = apfel::NF(std::sqrt(q2_max_val), thresholds);

    // Derive channels: 3 per quark (qq, gq, qg)
    auto channels = pineapfel::derive_channels(pineapfel::ProcessType::SIDIS,
        pineapfel::Observable::F2,
        pineapfel::Current::NC,
        pineapfel::CCSign::Plus,
        nf_max);

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
                    double as_power = std::pow(as_val, alpha_s);

                    const apfel::DoubleObject<apfel::Operator> *coeff = nullptr;
                    if (alpha_s == 0 && channel_type == 0) coeff = &sobj.C20qq;
                    else if (alpha_s == 1 && channel_type == 0)
                        coeff = &sobj.C21qq;
                    else if (alpha_s == 1 && channel_type == 1)
                        coeff = &sobj.C21gq;
                    else if (alpha_s == 1 && channel_type == 2)
                        coeff = &sobj.C21qg;

                    if (coeff == nullptr) continue;

                    const auto &terms = coeff->GetTerms();
                    for (const auto &term : terms) {
                        for (size_t ic = 0; ic < ch.pid_combinations.size();
                             ic++) {
                            int pid_pdf = ch.pid_combinations[ic][0];
                            int pid_ff  = ch.pid_combinations[ic][1];

                            apfel::Distribution pdf_dist(g,
                                [&](double const &xx) -> double {
                                    return toy_f(pid_pdf, xx);
                                });
                            apfel::Distribution ff_dist(g,
                                [&](double const &zz) -> double {
                                    return toy_f(pid_ff, zz);
                                });

                            double              val_x =
                                (term.object1 * pdf_dist).Evaluate(x_center);
                            double val_z =
                                (term.object2 * ff_dist).Evaluate(z_center);

                            result[ibin] += as_power * e_q_sq * ch.factors[ic] *
                                            term.coefficient * val_x * val_z;
                        }
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
    auto   sobj       = pineapfel::init_sidis_pol(g, thresholds);
    size_t nbins      = bin_x_bounds.size();

    // Determine nf_max
    double q2_max_val = 0;
    for (double q2 : q2_nodes) q2_max_val = std::max(q2_max_val, q2);
    int  nf_max   = apfel::NF(std::sqrt(q2_max_val), thresholds);

    // Derive channels: 3 per quark (qq, gq, qg)
    auto channels = pineapfel::derive_channels(pineapfel::ProcessType::SIDIS,
        pineapfel::Observable::F2,
        pineapfel::Current::NC,
        pineapfel::CCSign::Plus,
        nf_max);

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
                    double as_power = std::pow(as_val, alpha_s);

                    const apfel::DoubleObject<apfel::Operator> *coeff = nullptr;
                    if (alpha_s == 0 && channel_type == 0) coeff = &sobj.G10qq;
                    else if (alpha_s == 1 && channel_type == 0)
                        coeff = &sobj.G11qq;
                    else if (alpha_s == 1 && channel_type == 1)
                        coeff = &sobj.G11gq;
                    else if (alpha_s == 1 && channel_type == 2)
                        coeff = &sobj.G11qg;

                    if (coeff == nullptr) continue;

                    const auto &terms = coeff->GetTerms();
                    for (const auto &term : terms) {
                        for (size_t ic = 0; ic < ch.pid_combinations.size();
                             ic++) {
                            int pid_pdf = ch.pid_combinations[ic][0];
                            int pid_ff  = ch.pid_combinations[ic][1];

                            apfel::Distribution pdf_dist(g,
                                [&](double const &xx) -> double {
                                    return toy_f(pid_pdf, xx);
                                });
                            apfel::Distribution ff_dist(g,
                                [&](double const &zz) -> double {
                                    return toy_f(pid_ff, zz);
                                });

                            double              val_x =
                                (term.object1 * pdf_dist).Evaluate(x_center);
                            double val_z =
                                (term.object2 * ff_dist).Evaluate(z_center);

                            result[ibin] += as_power * e_q_sq * ch.factors[ic] *
                                            term.coefficient * val_x * val_z;
                        }
                    }
                }
            }
        }
    }

    return result;
}
