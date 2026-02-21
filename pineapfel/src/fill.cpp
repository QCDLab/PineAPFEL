#include <apfel/apfelxx.h>
#include <fill.h>
#include <sidis_api.h>

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
    bool                                   polarized,
    const apfel::Grid                     &g,
    const std::vector<double>             &thresholds)
    -> std::function<apfel::StructureFunctionObjects(double const &,
        std::vector<double> const &)> {
    if (process == ProcessType::DIS) {
        if (polarized) {
            if (current == Current::CC)
                throw std::runtime_error("build_grid: CC current is not "
                                         "supported for polarized DIS");
            switch (observable) {
            case Observable::F2:
                return apfel::Initializeg1NCObjectsZM(g, thresholds);
            case Observable::FL:
                return apfel::InitializegLNCObjectsZM(g, thresholds);
            case Observable::F3:
                return apfel::Initializeg4NCObjectsZM(g, thresholds);
            }
        } else if (current == Current::NC) {
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
        if (polarized)
            throw std::runtime_error(
                "build_grid: polarized SIA is not supported");
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

// Select the appropriate DoubleObject<Operator> from SidisObjects
// for a given perturbative order, channel type (qq/gq/qg), and observable.
static const apfel::DoubleObject<apfel::Operator> *select_sidis_coeff(
    const pineapfel::SidisCoeffs &sobj,
    int                           alpha_s,
    int                           channel_type, // 0=qq, 1=gq, 2=qg
    Observable                    observable,
    int                           nf) {
    if (observable == Observable::F2) {
        if (alpha_s == 0) {
            if (channel_type == 0) return &sobj.C20qq;
            return nullptr; // gq, qg start at NLO
        } else if (alpha_s == 1) {
            if (channel_type == 0) return &sobj.C21qq;
            if (channel_type == 1) return &sobj.C21gq;
            if (channel_type == 2) return &sobj.C21qg;
        } else if (alpha_s == 2) {
            if (channel_type == 0) {
                auto it = sobj.C22qq.find(nf);
                if (it != sobj.C22qq.end()) return &it->second;
            }
            return nullptr; // gq, qg not available at NNLO
        }
    } else if (observable == Observable::FL) {
        if (alpha_s == 0) return nullptr; // No LO FL
        if (alpha_s == 1) {
            if (channel_type == 0) return &sobj.CL1qq;
            if (channel_type == 1) return &sobj.CL1gq;
            if (channel_type == 2) return &sobj.CL1qg;
        }
        // No NNLO FL in APFEL++
    }
    return nullptr;
}

// Select the appropriate polarized DoubleObject from SidisPolCoeffs.
static const apfel::DoubleObject<apfel::Operator> *select_sidis_pol_coeff(
    const pineapfel::SidisPolCoeffs &sobj,
    int                              alpha_s,
    int                              channel_type, // 0=qq, 1=gq, 2=qg
    Observable                       observable,
    int                              nf) {
    // Only F2 (G1) is available for polarized SIDIS in APFEL++
    if (observable != Observable::F2) return nullptr;
    if (alpha_s == 0) {
        if (channel_type == 0) return &sobj.G10qq;
        return nullptr; // gq, qg start at NLO
    } else if (alpha_s == 1) {
        if (channel_type == 0) return &sobj.G11qq;
        if (channel_type == 1) return &sobj.G11gq;
        if (channel_type == 2) return &sobj.G11qg;
    } else if (alpha_s == 2) {
        if (channel_type == 0) {
            auto it = sobj.G12qq.find(nf);
            if (it != sobj.G12qq.end()) return &it->second;
        }
        return nullptr; // gq, qg not available at NNLO
    }
    return nullptr;
}

static pineappl_grid *build_grid_sidis(const GridDef &grid_def_in,
    const TheoryCard                                 &theory,
    const OperatorCard                               &op_card) {
    if (grid_def_in.observable == Observable::F3)
        throw std::runtime_error(
            "build_grid_sidis: F3 is not supported for SIDIS");
    if (grid_def_in.current == Current::CC)
        throw std::runtime_error(
            "build_grid_sidis: CC current is not supported for SIDIS");
    if (grid_def_in.polarized && grid_def_in.observable == Observable::FL)
        throw std::runtime_error(
            "build_grid_sidis: FL is not supported for polarized SIDIS");

    std::cout << "Building APFEL++ "
              << (grid_def_in.polarized ? "polarized " : "")
              << "SIDIS coefficient function grid..." << std::endl;

    // Auto-derive channels
    GridDef grid_def = grid_def_in;
    double  q2_max   = 0;
    for (const auto &bin : grid_def.bins)
        q2_max = std::max(q2_max, bin.upper[0]);
    int nf_max        = apfel::NF(std::sqrt(q2_max), theory.quark_thresholds);
    grid_def.channels = derive_channels(grid_def.process,
        grid_def.observable,
        grid_def.current,
        grid_def.cc_sign,
        nf_max);
    std::cout << "  Auto-derived " << grid_def.channels.size()
              << " channels for nf_max=" << nf_max << std::endl;

    // 1. Build APFEL++ x-space grid (same grid for x and z)
    std::vector<apfel::SubGrid> subgrids;
    for (const auto &sg : op_card.xgrid)
        subgrids.emplace_back(sg.n_knots, sg.x_min, sg.poly_degree);
    const apfel::Grid g{subgrids};

    // 2. Initialize SIDIS coefficient functions and build a unified selector.
    // The selector lambda captures the right coefficient struct by value so
    // that the fill loop below is identical for polarized and unpolarized
    // cases.
    using CoeffPtr = const apfel::DoubleObject<apfel::Operator> *;
    std::function<CoeffPtr(int, int, Observable, int)> get_coeff;
    if (grid_def_in.polarized) {
        auto sobj = std::make_shared<SidisPolCoeffs>(
            init_sidis_pol(g, theory.quark_thresholds));
        get_coeff = [sobj](int     alpha_s,
                        int        channel_type,
                        Observable obs,
                        int        nf) -> CoeffPtr {
            return select_sidis_pol_coeff(*sobj,
                alpha_s,
                channel_type,
                obs,
                nf);
        };
    } else {
        auto sobj = std::make_shared<SidisCoeffs>(
            init_sidis(g, theory.quark_thresholds));
        get_coeff = [sobj](int     alpha_s,
                        int        channel_type,
                        Observable obs,
                        int        nf) -> CoeffPtr {
            return select_sidis_coeff(*sobj, alpha_s, channel_type, obs, nf);
        };
    }

    // 3. Create empty PineAPPL grid
    pineappl_grid      *grid           = create_grid(grid_def);

    // 4. Determine grid nodes
    const auto         &joint_grid_vec = g.GetJointGrid().GetGrid();
    std::vector<double> x_nodes(joint_grid_vec.begin(), joint_grid_vec.end());
    const std::size_t   nx = x_nodes.size();
    // Same grid for z
    const std::size_t   nz = nx;

    std::vector<double> q2_nodes =
        derive_q2_nodes(grid_def.bins, theory.quark_thresholds);
    const std::size_t nq = q2_nodes.size();

    std::cout << "  Grid nodes: " << nq << " Q^2 x " << nx << " x x " << nz
              << " z points" << std::endl;

    // node_values: [q2_0..q2_{nq-1}, x_0..x_{nx-1}, z_0..z_{nz-1}]
    std::vector<double> node_values;
    node_values.reserve(nq + nx + nz);
    node_values.insert(node_values.end(), q2_nodes.begin(), q2_nodes.end());
    node_values.insert(node_values.end(), x_nodes.begin(), x_nodes.end());
    node_values.insert(node_values.end(), x_nodes.begin(), x_nodes.end());

    std::vector<std::size_t> shape = {nq, nx, nz};

    // 5. Precompute electroweak charges per Q² node
    struct Q2DataSidis {
        int                 nf;
        std::vector<double> charges;
    };
    std::vector<Q2DataSidis> q2_data(nq);
    for (std::size_t iq = 0; iq < nq; iq++) {
        double Q       = std::sqrt(q2_nodes[iq]);
        q2_data[iq].nf = apfel::NF(Q, theory.quark_thresholds);
        q2_data[iq].charges =
            apfel::ElectroWeakCharges(Q, false /* spacelike */);
    }

    // 6. Fill subgrids for each (order, channel, bin)
    // Channels are grouped: for quark q, indices are 3*(q-1)+type
    //   type: 0=qq, 1=gq, 2=qg
    for (std::size_t iord = 0; iord < grid_def.orders.size(); iord++) {
        int alpha_s = grid_def.orders[iord].alpha_s;
        if (alpha_s > 2) {
            std::cerr << "  Warning: skipping order alpha_s^" << alpha_s
                      << " (beyond NNLO)" << std::endl;
            continue;
        }

        for (std::size_t ich = 0; ich < grid_def.channels.size(); ich++) {
            int quark_idx    = static_cast<int>(ich / 3); // 0-based quark index
            int channel_type = static_cast<int>(ich % 3); // 0=qq, 1=gq, 2=qg

            for (std::size_t ibin = 0; ibin < grid_def.bins.size(); ibin++) {
                // Bin centers: dim 0=Q², dim 1=x, dim 2=z
                double              x_lo     = grid_def.bins[ibin].lower[1];
                double              x_hi     = grid_def.bins[ibin].upper[1];
                double              x_center = std::sqrt(x_lo * x_hi);

                double              z_lo     = grid_def.bins[ibin].lower[2];
                double              z_hi     = grid_def.bins[ibin].upper[2];
                double              z_center = std::sqrt(z_lo * z_hi);

                std::vector<double> subgrid(nq * nx * nz, 0.0);

                for (std::size_t iq = 0; iq < nq; iq++) {
                    int nf = q2_data[iq].nf;

                    // Skip if quark is not active at this Q²
                    if (quark_idx + 1 > nf) continue;

                    const auto *coeff = get_coeff(alpha_s,
                        channel_type,
                        grid_def.observable,
                        nf);
                    if (coeff == nullptr) continue;

                    // e_q² weight for this quark
                    double      e_q_sq = q2_data[iq].charges[quark_idx];

                    // Evaluate via outer product of the DoubleObject terms
                    const auto &terms  = coeff->GetTerms();
                    for (const auto &term : terms) {
                        double              c = term.coefficient;
                        apfel::Distribution dist_x =
                            term.object1.Evaluate(x_center);
                        apfel::Distribution dist_z =
                            term.object2.Evaluate(z_center);

                        const auto &vx = dist_x.GetDistributionJointGrid();
                        const auto &vz = dist_z.GetDistributionJointGrid();

                        for (std::size_t ix = 0; ix < nx && ix < vx.size();
                             ix++) {
                            for (std::size_t iz = 0; iz < nz && iz < vz.size();
                                 iz++) {
                                subgrid[iq * nx * nz + ix * nz + iz] +=
                                    e_q_sq * c * vx[ix] * vz[iz];
                            }
                        }
                    }
                }

                set_subgrid(grid, ibin, iord, ich, node_values, subgrid, shape);
            }
        }
    }

    std::cout << "SIDIS grid filled successfully." << std::endl;
    return grid;
}

pineappl_grid *build_grid(const GridDef &grid_def_in,
    const TheoryCard                    &theory,
    const OperatorCard                  &op_card) {
    if (grid_def_in.process == ProcessType::SIDIS)
        return build_grid_sidis(grid_def_in, theory, op_card);

    std::cout << "Building APFEL++ coefficient function grid..." << std::endl;

    // Auto-derive channels from observable and number of active flavors
    GridDef grid_def = grid_def_in;
    double  q2_max   = 0;
    for (const auto &bin : grid_def.bins)
        q2_max = std::max(q2_max, bin.upper[0]);
    int nf_max        = apfel::NF(std::sqrt(q2_max), theory.quark_thresholds);
    grid_def.channels = derive_channels(grid_def.process,
        grid_def.observable,
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
        grid_def.polarized,
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
