#include <evolution.h>
#include <apfel/apfelxx.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace pineapfel {

// ---- Internal helpers (not exposed) ----

static std::vector<std::size_t> unravel_index(std::size_t flat_index, const std::vector<std::size_t>& shape) {
    std::size_t ndim = shape.size();
    std::vector<std::size_t> coords(ndim);
    for (int i = ndim - 1; i >= 0; --i) {
        coords[i] = flat_index % shape[i];
        flat_index /= shape[i];
    }
    return coords;
}

// ---- QCD internal params and callback ----

struct ApfelxxParamsQCD {
    double mu0;
    int pert_ord;
    std::vector<pineappl_conv_type> conv_types;
    std::vector<double> thresholds;
    std::function<double(double const&)> alphas;
    std::map<pineappl_conv_type, std::map<int, apfel::DglapObjects>> DglapObj;
};

extern "C" void evolution_operators_qcd(
    std::size_t op_index,
    double fac1,
    const int* pids_in,
    const double* x_in,
    const int* pids_out,
    const double* x_out,
    const std::size_t* eko_shape,
    double* eko_buffer,
    void* params
) {
    ApfelxxParamsQCD* op_params = static_cast<ApfelxxParamsQCD*>(params);
    pineappl_conv_type conv_type = op_params->conv_types[op_index];
    std::cout << "Convolution type: " << conv_type << ", mu = " << sqrt(fac1) << " GeV" << std::endl;

    std::vector<std::size_t> shape(eko_shape, eko_shape + 4);
    std::size_t flat_len = std::accumulate(shape.begin(), shape.end(), std::size_t{1}, std::multiplies<std::size_t>());

    // Zero operator
    const apfel::Operator ZeroOp{
        op_params->DglapObj.at(conv_type).begin()->second.SplittingFunctions.begin()->second.GetObjects().begin()->second.GetGrid(),
        apfel::Null{}, apfel::eps5
    };

    // Build evolution operator and rotate to physical basis (13x13)
    std::map<std::pair<int, int>, apfel::Operator> EvOpPhys;
    const std::map<int, apfel::Operator> EvOpPlusMinus =
        BuildDglap(op_params->DglapObj.at(conv_type), op_params->mu0, op_params->pert_ord, op_params->alphas)
            ->Evaluate(sqrt(fac1)).GetObjects();

    for (int i = 0; i < 13; i++) {
        for (int j = 0; j < 13; j++) {
            EvOpPhys.insert({{i, j}, ZeroOp});
            for (int k = 0; k < 13; k++) {
                for (int l = 0; l < 13; l++) {
                    if (apfel::GkjPhys.count({k, l}) == 0) continue;
                    EvOpPhys.at({i, j}) += (apfel::RotPlusMinusToPhys[i][k] * apfel::RotPhysToPlusMinus[l][j])
                                           * EvOpPlusMinus.at(apfel::GkjPhys.at({k, l}));
                }
            }
        }
    }

    // Flatten into eko_buffer
    for (std::size_t i = 0; i != flat_len; i++) {
        std::vector<std::size_t> coords = unravel_index(i, shape);

        // Skip photon
        if (pids_in[coords[0]] == 22 || pids_out[coords[2]] == 22) continue;

        const std::pair<int, int> ij{
            (pids_in[coords[0]] == 21 ? 0 : pids_in[coords[0]]) + 6,
            (pids_out[coords[2]] == 21 ? 0 : pids_out[coords[2]]) + 6
        };
        eko_buffer[i] = EvOpPhys.at(ij).Evaluate(x_in[coords[1]]).GetDistributionJointGrid()[coords[3]]
                        * x_out[coords[3]] / x_in[coords[1]];
    }
}

// ---- QED internal params and callback ----

struct ApfelxxParamsQED {
    double mu0;
    int pert_ord;
    std::vector<pineappl_conv_type> conv_types;
    std::vector<double> quark_thresholds;
    std::vector<double> lepton_thresholds;
    std::function<double(double const&)> alphas;
    std::function<double(double const&)> alphaem;
    std::map<pineappl_conv_type, std::map<int, apfel::DglapObjectsQCDQED>> DglapObj;
};

extern "C" void evolution_operators_qed(
    std::size_t op_index,
    double fac1,
    const int* pids_in,
    const double* x_in,
    const int* pids_out,
    const double* x_out,
    const std::size_t* eko_shape,
    double* eko_buffer,
    void* params
) {
    ApfelxxParamsQED* op_params = static_cast<ApfelxxParamsQED*>(params);
    pineappl_conv_type conv_type = op_params->conv_types[op_index];
    std::cout << "Convolution type: " << conv_type << ", mu = " << sqrt(fac1) << " GeV" << std::endl;

    std::vector<std::size_t> shape(eko_shape, eko_shape + 4);
    std::size_t flat_len = std::accumulate(shape.begin(), shape.end(), std::size_t{1}, std::multiplies<std::size_t>());

    // Zero operator
    const apfel::Operator ZeroOp{
        op_params->DglapObj.at(conv_type).begin()->second.SplittingFunctions.begin()->second.GetObjects().begin()->second.GetGrid(),
        apfel::Null{}, apfel::eps5
    };

    // Build evolution operator and rotate to physical basis (20x20)
    std::map<std::pair<int, int>, apfel::Operator> EvOpPhys;
    const std::map<int, apfel::Operator> EvOpPlusMinus =
        BuildDglap(op_params->DglapObj.at(conv_type), op_params->mu0, op_params->pert_ord, op_params->alphas, op_params->alphaem)
            ->Evaluate(sqrt(fac1)).GetObjects();

    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {
            EvOpPhys.insert({{i, j}, ZeroOp});
            for (int k = 0; k < 20; k++) {
                for (int l = 0; l < 20; l++) {
                    if (apfel::GkjQCDQED.count({k, l}) == 0) continue;
                    EvOpPhys.at({i, j}) += (apfel::RotPlusMinusQCDQEDToPhys[i][k] * apfel::RotPhysToPlusMinusQCDQED[l][j])
                                           * EvOpPlusMinus.at(apfel::GkjQCDQED.at({k, l}));
                }
            }
        }
    }

    // Flatten into eko_buffer
    for (std::size_t i = 0; i != flat_len; i++) {
        std::vector<std::size_t> coords = unravel_index(i, shape);

        // Skip leptons
        if (std::abs(pids_in[coords[0]]) == 11 || std::abs(pids_in[coords[0]]) == 13 || std::abs(pids_in[coords[0]]) == 15 ||
            std::abs(pids_out[coords[2]]) == 11 || std::abs(pids_out[coords[2]]) == 13 || std::abs(pids_out[coords[2]]) == 15)
            continue;

        int ii = -1;
        if (pids_in[coords[0]] == 21)
            ii = 6;
        else if (pids_in[coords[0]] == 22)
            ii = 16;
        else
            ii = pids_in[coords[0]] + 6;

        int jj = -1;
        if (pids_out[coords[2]] == 21)
            jj = 6;
        else if (pids_out[coords[2]] == 22)
            jj = 16;
        else
            jj = pids_out[coords[2]] + 6;

        eko_buffer[i] = EvOpPhys.at({ii, jj}).Evaluate(x_in[coords[1]]).GetDistributionJointGrid()[coords[3]]
                        * x_out[coords[3]] / x_in[coords[1]];
    }
}

// ---- Public evolve function ----

pineappl_grid* evolve(
    pineappl_grid* grid,
    const TheoryCard& theory,
    const OperatorCard& op_card
) {
    // Verify PID basis
    const pineappl_pid_basis pid_basis = pineappl_grid_pid_basis(grid);
    assert(pid_basis == PINEAPPL_PID_BASIS_PDG);

    // Get convolution types from the grid
    std::size_t n_convs = pineappl_grid_convolutions_len(grid);
    std::vector<pineappl_conv_type> conv_types(n_convs);
    pineappl_grid_conv_types(grid, conv_types.data());

    // Unique convolution types
    std::vector<pineappl_conv_type> unique_convs;
    for (std::size_t i = 0; i != n_convs; i++) {
        pineappl_conv_type conv = conv_types[i];
        if (std::find(unique_convs.begin(), unique_convs.end(), conv) == unique_convs.end()) {
            unique_convs.push_back(conv);
        }
    }

    // Get evolve info shape and values
    std::vector<std::size_t> evinfo_shape(5);
    pineappl_grid_evolve_info_shape(grid, nullptr, evinfo_shape.data());

    std::vector<double> fac1(evinfo_shape[0]);
    std::vector<double> frg1(evinfo_shape[1]);
    std::vector<int> pids_in(evinfo_shape[2]);
    std::vector<double> x_in(evinfo_shape[3]);
    std::vector<double> ren1(evinfo_shape[4]);
    pineappl_grid_evolve_info(grid, nullptr, fac1.data(), frg1.data(), pids_in.data(), x_in.data(), ren1.data());

    // Build APFEL++ x-space grid
    std::vector<apfel::SubGrid> subgrids;
    for (const auto& sg : op_card.xgrid) {
        subgrids.emplace_back(sg.n_knots, sg.x_min, sg.poly_degree);
    }
    const apfel::Grid g{subgrids};

    // Get joint grid for x_out
    std::vector<double> x_out(g.GetJointGrid().nx() + 1);
    for (std::size_t i = 0; i < x_out.size(); i++)
        x_out[i] = g.GetJointGrid().GetGrid()[i];

    // Output PIDs from theory card
    std::vector<int> pids_out = theory.flavors;

    // EKO shape
    const std::vector<std::size_t> shape = {pids_in.size(), x_in.size(), pids_out.size(), x_out.size()};

    // Construct operator info slices
    std::vector<pineappl_conv_type> convtypes(unique_convs.size());
    std::vector<pineappl_operator_info> opinfo_slices(unique_convs.size() * fac1.size());
    for (std::size_t i = 0; i != unique_convs.size(); i++) {
        for (std::size_t j = 0; j != fac1.size(); j++) {
            pineappl_operator_info opinfo = {
                std::pow(theory.mu0, 2),
                fac1[j],
                pid_basis,
                unique_convs[i],
            };
            opinfo_slices[i * fac1.size() + j] = opinfo;
        }
        convtypes[i] = unique_convs[i];
    }

    // Tabulation params
    const auto& tab = op_card.tabulation;

    pineappl_grid* fktable = nullptr;

    if (theory.qed) {
        // ---- QED path ----
        auto* op_params = new ApfelxxParamsQED;
        op_params->mu0 = theory.mu0;
        op_params->pert_ord = theory.pert_ord;
        op_params->quark_thresholds = theory.quark_thresholds;
        op_params->lepton_thresholds = theory.lepton_thresholds;
        op_params->conv_types = convtypes;

        // Construct coupled alpha_s + alpha_em
        apfel::AlphaQCDQED a{theory.alpha_qcd_ref, theory.alpha_qed_ref, theory.q_ref,
                             theory.quark_thresholds, theory.lepton_thresholds, theory.pert_ord};
        const apfel::TabulateObject<apfel::matrix<double>> Couplings{a, tab.n_points, tab.q_min, static_cast<double>(tab.n_steps), tab.interp_degree};
        op_params->alphas  = [=](double const& mu) -> double { return Couplings.Evaluate(mu)(0, 0); };
        op_params->alphaem = [=](double const& mu) -> double { return Couplings.Evaluate(mu)(1, 0); };

        // Initialize DGLAP objects for QCD+QED (only UNPOL_PDF supported for QED)
        op_params->DglapObj.insert({
            pineappl_conv_type::PINEAPPL_CONV_TYPE_UNPOL_PDF,
            InitializeDglapObjectsQCDQED(g, theory.quark_thresholds, theory.lepton_thresholds, true)
        });

        // Alpha_s table for renormalization scales
        std::vector<double> alphas_table;
        for (double q2 : ren1) {
            alphas_table.push_back(op_params->alphas(sqrt(q2)));
        }

        fktable = pineappl_grid_evolve(
            grid,
            unique_convs.size(),
            evolution_operators_qed,
            opinfo_slices.data(),
            pids_in.data(),
            x_in.data(),
            pids_out.data(),
            x_out.data(),
            shape.data(),
            static_cast<void*>(op_params),
            nullptr,
            op_card.xi.data(),
            ren1.data(),
            alphas_table.data()
        );

        delete op_params;
    } else {
        // ---- QCD path ----
        auto* op_params = new ApfelxxParamsQCD;
        op_params->mu0 = theory.mu0;
        op_params->pert_ord = theory.pert_ord;
        op_params->thresholds = theory.quark_thresholds;
        op_params->conv_types = convtypes;

        // Construct running strong coupling
        apfel::AlphaQCD a{theory.alpha_qcd_ref, theory.q_ref, theory.quark_thresholds, theory.pert_ord};
        const apfel::TabulateObject<double> Alphas{a, tab.n_points, tab.q_min, static_cast<double>(tab.n_steps), tab.interp_degree};
        op_params->alphas = [=](double const& mu) -> double { return Alphas.Evaluate(mu); };

        // Initialize all possible DGLAP objects
        op_params->DglapObj.insert({
            pineappl_conv_type::PINEAPPL_CONV_TYPE_UNPOL_PDF,
            InitializeDglapObjectsQCDPhys(g, theory.quark_thresholds, true)
        });
        op_params->DglapObj.insert({
            pineappl_conv_type::PINEAPPL_CONV_TYPE_POL_PDF,
            InitializeDglapObjectsQCDpolPhys(g, theory.quark_thresholds, true)
        });
        op_params->DglapObj.insert({
            pineappl_conv_type::PINEAPPL_CONV_TYPE_UNPOL_FF,
            InitializeDglapObjectsQCDTPhys(g, theory.quark_thresholds, true)
        });

        // Alpha_s table for renormalization scales
        std::vector<double> alphas_table;
        for (double q2 : ren1) {
            alphas_table.push_back(op_params->alphas(sqrt(q2)));
        }

        fktable = pineappl_grid_evolve(
            grid,
            unique_convs.size(),
            evolution_operators_qcd,
            opinfo_slices.data(),
            pids_in.data(),
            x_in.data(),
            pids_out.data(),
            x_out.data(),
            shape.data(),
            static_cast<void*>(op_params),
            nullptr,
            op_card.xi.data(),
            ren1.data(),
            alphas_table.data()
        );

        delete op_params;
    }

    return fktable;
}

} // namespace pineapfel
