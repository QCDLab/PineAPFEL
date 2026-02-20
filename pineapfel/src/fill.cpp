#include <apfel/apfelxx.h>
#include <fill.h>

#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

namespace pineapfel {

// Select the appropriate APFEL++ structure function initializer.
static auto select_initializer(ProcessType process,
    Observable                             observable,
    Current                                current,
    CCSign                                 cc_sign,
    const apfel::Grid                     &g,
    const std::vector<double>             &thresholds)
    -> std::function<apfel::StructureFunctionObjects(double const &,
        std::vector<double> const &)> {
    if (process == ProcessType::DIS) {
        if (current == Current::NC) {
            switch (observable) {
            case Observable::F2:
                return apfel::InitializeF2NCObjectsZM(g, thresholds);
            case Observable::FL:
                return apfel::InitializeFLNCObjectsZM(g, thresholds);
            case Observable::F3:
                return apfel::InitializeF3NCObjectsZM(g, thresholds);
            }
        } else {
            // CC
            if (cc_sign == CCSign::Plus) {
                switch (observable) {
                case Observable::F2:
                    return apfel::InitializeF2CCPlusObjectsZM(g, thresholds);
                case Observable::FL:
                    return apfel::InitializeFLCCPlusObjectsZM(g, thresholds);
                case Observable::F3:
                    return apfel::InitializeF3CCPlusObjectsZM(g, thresholds);
                }
            } else {
                switch (observable) {
                case Observable::F2:
                    return apfel::InitializeF2CCMinusObjectsZM(g, thresholds);
                case Observable::FL:
                    return apfel::InitializeFLCCMinusObjectsZM(g, thresholds);
                case Observable::F3:
                    return apfel::InitializeF3CCMinusObjectsZM(g, thresholds);
                }
            }
        }
    } else if (process == ProcessType::SIA) {
        if (current == Current::CC)
            throw std::runtime_error(
                "build_grid: CC current is not supported for SIA");
        switch (observable) {
        case Observable::F2:
            return apfel::InitializeF2NCObjectsZMT(g, thresholds);
        case Observable::FL:
            return apfel::InitializeFLNCObjectsZMT(g, thresholds);
        case Observable::F3:
            return apfel::InitializeF3NCObjectsZMT(g, thresholds);
        }
    }

    throw std::runtime_error(
        "build_grid: unsupported process/observable combination");
}

// Collect unique Q^2 nodes from bin edges with geometric intermediate points.
static std::vector<double> derive_q2_nodes(const std::vector<BinDef> &bins,
    const std::vector<double> &thresholds,
    int                        n_intermediate = 3) {
    std::set<double> q2_set;

    for (const auto &bin : bins) {
        double q2_lo = bin.lower[0];
        double q2_hi = bin.upper[0];
        q2_set.insert(q2_lo);
        q2_set.insert(q2_hi);

        if (n_intermediate > 0) {
            double log_lo = std::log(q2_lo);
            double log_hi = std::log(q2_hi);
            for (int i = 1; i <= n_intermediate; i++) {
                double frac = static_cast<double>(i) / (n_intermediate + 1);
                q2_set.insert(std::exp(log_lo + frac * (log_hi - log_lo)));
            }
        }
    }

    double q2_min = *q2_set.begin();
    double q2_max = *q2_set.rbegin();
    for (double thr : thresholds) {
        double q2_thr = thr * thr;
        if (q2_thr > q2_min && q2_thr < q2_max) q2_set.insert(q2_thr);
    }

    return std::vector<double>(q2_set.begin(), q2_set.end());
}

// Build the coefficient function operator for a given channel.
//
// Physical-basis decomposition:
//   F = sum_q C_q (x) (q +/- qbar) + C_g (x) g
// with:
//   C_q = w_q * CNS + (SumW / 6) * (CS - CNS)
//   C_g = SumW * CG
// where w_q is the per-quark weight (electroweak charge for NC, CKM weight
// for CC), SumW = sum of w_q for nf active flavors, and the factor 6 matches
// the hardcoded normalization in APFEL++'s DISNCBasis/DISCCBasis.
// APFEL++ sets CS = CNS and CG = 0 where the physics requires it, so this
// general formula works for all observables and currents.
static apfel::Operator build_channel_operator(const ChannelDef &channel,
    const std::map<int, apfel::Operator>                       &ops,
    const std::vector<double>                                  &charges,
    int                                                         nf) {
    const apfel::Operator &CNS    = ops.at(apfel::DISNCBasis::CNS);
    const apfel::Operator &CS     = ops.at(apfel::DISNCBasis::CS);
    const apfel::Operator &CG     = ops.at(apfel::DISNCBasis::CG);

    double                 sum_ch = 0;
    for (int i = 0; i < nf && i < static_cast<int>(charges.size()); i++)
        sum_ch += charges[i];

    // Detect gluon channel: any PID combination containing only PID 21
    bool is_gluon = false;
    for (const auto &combo : channel.pid_combinations) {
        if (combo.size() == 1 && combo[0] == 21) {
            is_gluon = true;
            break;
        }
    }

    if (is_gluon) { return sum_ch * CG; }

    // Quark channel: extract the quark PID
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

    // C_q = w_q * CNS + (sum_w / 6) * (CS - CNS)
    // Works for both NC and CC: APFEL++ sets operators to zero where needed
    return e_q_sq * CNS + (sum_ch / 6.0) * (CS - CNS);
}

pineappl_grid *build_grid(const GridDef &grid_def_in,
    const TheoryCard                    &theory,
    const OperatorCard                  &op_card) {
    if (grid_def_in.process == ProcessType::SIDIS)
        throw std::runtime_error("build_grid: SIDIS is not supported");

    std::cout << "Building APFEL++ coefficient function grid..." << std::endl;

    // Auto-derive channels from observable and number of active flavors
    GridDef grid_def = grid_def_in;
    double  q2_max   = 0;
    for (const auto &bin : grid_def.bins)
        q2_max = std::max(q2_max, bin.upper[0]);
    int nf_max        = apfel::NF(std::sqrt(q2_max), theory.quark_thresholds);
    grid_def.channels = derive_channels(grid_def.observable,
        grid_def.current,
        grid_def.cc_sign,
        nf_max);
    std::cout << "  Auto-derived " << grid_def.channels.size()
              << " channels for nf_max=" << nf_max << std::endl;

    // 1. Build APFEL++ x-space grid
    std::vector<apfel::SubGrid> subgrids;
    for (const auto &sg : op_card.xgrid)
        subgrids.emplace_back(sg.n_knots, sg.x_min, sg.poly_degree);
    const apfel::Grid   g{subgrids};

    // 2. Initialize structure function objects
    auto                sf_init        = select_initializer(grid_def.process,
        grid_def.observable,
        grid_def.current,
        grid_def.cc_sign,
        g,
        theory.quark_thresholds);

    // 3. Create empty PineAPPL grid
    pineappl_grid      *grid           = create_grid(grid_def);

    // 4. Determine grid nodes
    const auto         &joint_grid_vec = g.GetJointGrid().GetGrid();
    std::vector<double> x_nodes(joint_grid_vec.begin(), joint_grid_vec.end());
    const std::size_t   nx = x_nodes.size();

    std::vector<double> q2_nodes =
        derive_q2_nodes(grid_def.bins, theory.quark_thresholds);
    const std::size_t nq = q2_nodes.size();

    std::cout << "  Grid nodes: " << nq << " Q^2 x " << nx << " x/z points"
              << std::endl;

    // Concatenated node_values: [q2_0..q2_{nq-1}, x_0..x_{nx-1}]
    std::vector<double> node_values;
    node_values.reserve(nq + nx);
    node_values.insert(node_values.end(), q2_nodes.begin(), q2_nodes.end());
    node_values.insert(node_values.end(), x_nodes.begin(), x_nodes.end());

    std::vector<std::size_t> shape    = {nq, nx};

    bool                     timelike = (grid_def.process == ProcessType::SIA);
    bool                     is_cc    = (grid_def.current == Current::CC);

    // 5. Precompute structure function objects and operators for each Q^2 node.
    //    Structure: sf_data[iq] = { order -> {CNS, CS, CG} operators }
    //    Also store nf and charges per Q^2 node.
    struct Q2Data {
        int                 nf;
        std::vector<double> charges;
        std::map<int, std::map<int, apfel::Operator>>
            order_ops; // order -> operator map
    };

    std::vector<Q2Data> q2_data(nq);

    for (std::size_t iq = 0; iq < nq; iq++) {
        double Q       = std::sqrt(q2_nodes[iq]);
        q2_data[iq].nf = apfel::NF(Q, theory.quark_thresholds);

        // For NC: use electroweak charges as per-quark weights
        // For CC: compute per-quark CKM weights and pass CKM vector to
        //         APFEL++ initializer
        std::vector<double> init_charges;
        if (is_cc) {
            // Per-quark CKM weight: sum of CKM^2 elements where quark
            // participates, filtered by active partner flavors.
            // Divided by 2 because CC Plus/Minus are defined as
            // (F(nu) +/- F(nubar)) / 2, and each CKM element contributes
            // to both the up-type and down-type quark channels.
            int                 nf = q2_data[iq].nf;
            std::vector<double> weights(nf, 0.0);
            for (int q = 1; q <= nf; q++) {
                bool is_down = (q % 2 == 1);
                if (is_down) {
                    int d_gen = (q + 1) / 2; // d->1, s->2, b->3
                    for (int u_gen = 1; u_gen <= 3; u_gen++) {
                        int partner_pid = 2 * u_gen; // u->2, c->4, t->6
                        if (partner_pid <= nf)
                            weights[q - 1] +=
                                theory.ckm[(u_gen - 1) * 3 + (d_gen - 1)];
                    }
                } else {
                    int u_gen = q / 2; // u->1, c->2, t->3
                    for (int d_gen = 1; d_gen <= 3; d_gen++) {
                        int partner_pid = 2 * d_gen - 1; // d->1, s->3, b->5
                        if (partner_pid <= nf)
                            weights[q - 1] +=
                                theory.ckm[(u_gen - 1) * 3 + (d_gen - 1)];
                    }
                }
                weights[q - 1] /= 2.0;
            }
            q2_data[iq].charges = weights;
            init_charges        = theory.ckm;
        } else {
            q2_data[iq].charges = apfel::ElectroWeakCharges(Q, timelike);
            init_charges        = q2_data[iq].charges;
        }

        auto FObjQ = sf_init(Q, init_charges);

        // Extract operators at each perturbative order using k=1
        // (operators are the same for any k; only ConvBasis differs)
        if (FObjQ.C0.count(1))
            q2_data[iq].order_ops[0] = FObjQ.C0.at(1).GetObjects();
        if (FObjQ.C1.count(1))
            q2_data[iq].order_ops[1] = FObjQ.C1.at(1).GetObjects();
        if (FObjQ.C2.count(1))
            q2_data[iq].order_ops[2] = FObjQ.C2.at(1).GetObjects();
    }

    // 6. Fill subgrids for each (bin, order, channel)
    for (std::size_t iord = 0; iord < grid_def.orders.size(); iord++) {
        int alpha_s = grid_def.orders[iord].alpha_s;
        if (alpha_s > 2) {
            std::cerr << "  Warning: skipping order alpha_s^" << alpha_s
                      << " (beyond NNLO)" << std::endl;
            continue;
        }

        for (std::size_t ich = 0; ich < grid_def.channels.size(); ich++) {
            for (std::size_t ibin = 0; ibin < grid_def.bins.size(); ibin++) {
                // Bin center in x/z dimension (last dimension, geometric mean)
                double              x_lo     = grid_def.bins[ibin].lower.back();
                double              x_hi     = grid_def.bins[ibin].upper.back();
                double              x_center = std::sqrt(x_lo * x_hi);

                // Build the full [nq, nx] subgrid
                std::vector<double> subgrid(nq * nx, 0.0);

                for (std::size_t iq = 0; iq < nq; iq++) {
                    if (q2_data[iq].order_ops.count(alpha_s) == 0) continue;

                    apfel::Operator C_channel =
                        build_channel_operator(grid_def.channels[ich],
                            q2_data[iq].order_ops[alpha_s],
                            q2_data[iq].charges,
                            q2_data[iq].nf);

                    // Evaluate operator at x_center -> Distribution on joint
                    // grid
                    apfel::Distribution dist = C_channel.Evaluate(x_center);
                    const std::vector<double> &vals =
                        dist.GetDistributionJointGrid();

                    for (std::size_t ix = 0; ix < nx && ix < vals.size(); ix++)
                        subgrid[iq * nx + ix] = vals[ix];
                }

                set_subgrid(grid, ibin, iord, ich, node_values, subgrid, shape);
            }
        }
    }

    std::cout << "Grid filled successfully." << std::endl;
    return grid;
}

} // namespace pineapfel
