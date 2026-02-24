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

using SFInitFn = std::function<apfel::StructureFunctionObjects(double const &,
    std::vector<double> const &)>;

struct WeightedSFInit {
    SFInitFn init;
    double   sign;  // +1 or -1
    int      actnf; // -1 = ZM (Q-dependent nf); >=0 = fixed (massive)
};

// Build the ZM (or polarized ZM) initializer — handles all
// process/observable/current combinations.
static SFInitFn make_zm_init(ProcessType process,
    Observable                           observable,
    Current                              current,
    CCSign                               cc_sign,
    bool                                 polarized,
    const apfel::Grid                   &g,
    const std::vector<double>           &thresholds) {
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

// Build the Massive or MassiveZero initializer for NC F2/FL.
// Returns nullptr for unsupported cases (F3, CC, polarized, SIA).
static SFInitFn make_massive_init(Observable obs,
    bool                                     zero,
    const apfel::Grid                       &g,
    const TheoryCard                        &theory) {
    const auto &M = theory.heavy_quark_masses;
    if (!zero) {
        if (obs == Observable::F2)
            return apfel::InitializeF2NCObjectsMassive(g,
                M,
                1e-5,
                theory.mass_nxi,
                theory.mass_ximin,
                theory.mass_ximax,
                theory.mass_intdeg,
                theory.mass_lambda,
                theory.mass_imod);
        else
            return apfel::InitializeFLNCObjectsMassive(g,
                M,
                1e-5,
                theory.mass_nxi,
                theory.mass_ximin,
                theory.mass_ximax,
                theory.mass_intdeg,
                theory.mass_lambda,
                theory.mass_imod);
    } else {
        if (obs == Observable::F2)
            return apfel::InitializeF2NCObjectsMassiveZero(g,
                M,
                1e-5,
                theory.mass_nxi,
                theory.mass_ximin,
                theory.mass_ximax,
                theory.mass_intdeg,
                theory.mass_lambda);
        else
            return apfel::InitializeFLNCObjectsMassiveZero(g,
                M,
                1e-5,
                theory.mass_nxi,
                theory.mass_ximin,
                theory.mass_ximax,
                theory.mass_intdeg,
                theory.mass_lambda);
    }
}

// Select the set of weighted initializers for the requested mass scheme.
// Returns a vector of {init, sign} pairs to be accumulated into operator maps.
static std::vector<WeightedSFInit> select_initializers(ProcessType process,
    Observable                                                     observable,
    Current                                                        current,
    CCSign                                                         cc_sign,
    bool                                                           polarized,
    MassScheme                                                     mass_scheme,
    const apfel::Grid                                             &g,
    const TheoryCard                                              &theory) {
    // Warn and fall back to ZM for cases where massive scheme is unsupported
    bool massive_unsupported =
        (current == Current::CC || polarized || process == ProcessType::SIA ||
            observable == Observable::F3);

    if (mass_scheme != MassScheme::ZM && massive_unsupported) {
        std::cerr << "Warning: FFN/FONLL not supported for this "
                     "process/observable/current; falling back to ZM.\n";
        mass_scheme = MassScheme::ZM;
    }

    auto zm            = make_zm_init(process,
        observable,
        current,
        cc_sign,
        polarized,
        g,
        theory.quark_thresholds);

    // Count light quarks (mass ≈ 0) to determine actnf for massive inits
    int  actnf_massive = 0;
    for (auto m : theory.heavy_quark_masses)
        if (m < 1e-8) actnf_massive++;

    switch (mass_scheme) {
    case MassScheme::ZM:
        return {
            {zm, +1.0, -1}
        };

    case MassScheme::FFN: {
        auto ffn = make_massive_init(observable, false, g, theory);
        return {
            {ffn, +1.0, actnf_massive}
        };
    }

    case MassScheme::MassiveZero:
        // APFEL++ InitializeF2NCObjectsMassiveZero sets the total channel (k=0)
        // to {CNS:Zero, CS:Zero, CG:Zero} ("Not") for all perturbative orders.
        // BuildStructureFunctions only uses k=0, so it gives BSF=0.
        // Return an empty initializer list so PineAPFEL also contributes zero.
        return {};

    case MassScheme::FONLL: {
        auto ffn = make_massive_init(observable, false, g, theory);
        // Standard FONLL = F_ZM + F_FFN - F_MassiveZero. APFEL++ sets the total
        // channel (k=0) of InitializeF2NCObjectsMassiveZero to zero ("Not") at
        // all orders, so F_MassiveZero = 0 here and the double-counting
        // subtraction is missing. The result is F_ZM + F_FFN, which over-counts
        // the light-quark ZM contribution.
        std::cerr
            << "Warning: FONLL is not fully implemented. "
               "The MassiveZero subtraction term (massless limit of F_FFN) "
               "is missing because APFEL++ does not expose it via the total "
               "channel (k=0). The grid will contain F_ZM + F_FFN instead "
               "of the correct F_ZM + F_FFN - F_MassiveZero.\n";
        return {
            { zm, +1.0,            -1},
            {ffn, +1.0, actnf_massive}
        };
    }
    }
    throw std::runtime_error("select_initializers: unhandled mass scheme");
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

// Select the appropriate DoubleObject<Operator> from SidisObjects for
// a given perturbative order, channel type (qq/gq/qg), and observable.
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

    // Infer polarization from the first convolution type (PDF slot):
    //   POL_PDF → polarized PDF coefficient functions
    //   POL_FF  → polarized FF (not currently supported, handled below)
    bool polarized =
        !grid_def_in.convolution_types.empty() &&
        (grid_def_in.convolution_types[0] == PINEAPPL_CONV_TYPE_POL_PDF ||
            grid_def_in.convolution_types[0] == PINEAPPL_CONV_TYPE_POL_FF);

    if (polarized && grid_def_in.observable == Observable::FL)
        throw std::runtime_error(
            "build_grid_sidis: FL is not supported for polarized SIDIS");

    std::cout << "Building APFEL++ " << (polarized ? "polarized " : "")
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
    if (polarized) {
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
    const apfel::Grid g{subgrids};

    // 2. Initialize structure function objects (one or more weighted inits).
    // Infer polarization from the first convolution type (PDF slot): POL_PDF or
    // POL_FF signals that polarized coefficient functions should be used.
    bool              polarized =
        !grid_def.convolution_types.empty() &&
        (grid_def.convolution_types[0] == PINEAPPL_CONV_TYPE_POL_PDF ||
            grid_def.convolution_types[0] == PINEAPPL_CONV_TYPE_POL_FF);
    auto                weighted_inits = select_initializers(grid_def.process,
        grid_def.observable,
        grid_def.current,
        grid_def.cc_sign,
        polarized,
        grid_def.mass_scheme,
        g,
        theory);

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

    std::vector<std::size_t> shape     = {nq, nx};

    bool                     timelike  = (grid_def.process == ProcessType::SIA);
    bool                     is_cc     = (grid_def.current == Current::CC);

    // Pre-determine gluon channel index (gluon is the last channel when
    // present)
    int                      gluon_ich = -1;
    bool                     has_gluon_channel  = false;
    int                      num_quark_channels = 0;
    for (int ich = 0; ich < (int)grid_def.channels.size(); ich++) {
        for (const auto &combo : grid_def.channels[ich].pid_combinations) {
            if (combo.size() == 1 && combo[0] == 21) {
                gluon_ich         = ich;
                has_gluon_channel = true;
                break;
            }
        }
        if (has_gluon_channel) break;
    }
    num_quark_channels =
        (int)grid_def.channels.size() - (has_gluon_channel ? 1 : 0);

    // 5. Precompute per-channel operators for each (Q^2, order) node.
    //    channel_ops[order][channel_idx] accumulates the combined operator
    //    from all weighted initializers.  Gluon and quark channels are
    //    handled separately to correctly incorporate heavy-quark
    //    contributions in the massive scheme.
    struct Q2Data {
        int                                           nf;
        // 6-entry EW charges (NC) or nf-entry CKM (CC):
        // channel_ops[order][channel_idx] = accumulated operator
        std::vector<double>                           charges;
        std::map<int, std::map<int, apfel::Operator>> channel_ops;
    };

    std::vector<Q2Data> q2_data(nq);

    // Helper: accumulate op*sign into target[ich].
    // NOTE: APFEL++ Operator::operator= is infinitely recursive (APFEL++ bug).
    // Use only emplace (first occurrence) and operator+= (subsequent) to avoid
    // it.
    auto add_to_channel = [](std::map<int, apfel::Operator> &target,
                              int                            ich,
                              const apfel::Operator         &op,
                              double                         sign) {
        auto it = target.find(ich);
        if (it == target.end()) target.emplace(ich, op * sign);
        else it->second += op * sign;
    };

    for (std::size_t iq = 0; iq < nq; iq++) {
        double Q       = std::sqrt(q2_nodes[iq]);
        q2_data[iq].nf = apfel::NF(Q, theory.quark_thresholds);

        // Compute per-quark charges and initializer charges
        std::vector<double> init_charges;
        if (is_cc) {
            int                 nf = q2_data[iq].nf;
            std::vector<double> weights(nf, 0.0);
            for (int q = 1; q <= nf; q++) {
                bool is_down = (q % 2 == 1);
                if (is_down) {
                    int d_gen = (q + 1) / 2;
                    for (int u_gen = 1; u_gen <= 3; u_gen++) {
                        int partner_pid = 2 * u_gen;
                        if (partner_pid <= nf)
                            weights[q - 1] +=
                                theory.ckm[(u_gen - 1) * 3 + (d_gen - 1)];
                    }
                } else {
                    int u_gen = q / 2;
                    for (int d_gen = 1; d_gen <= 3; d_gen++) {
                        int partner_pid = 2 * d_gen - 1;
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
            // ElectroWeakCharges always returns 6 entries
            q2_data[iq].charges = apfel::ElectroWeakCharges(Q, timelike);
            init_charges        = q2_data[iq].charges;
        }

        const std::vector<double> &charges = q2_data[iq].charges;

        for (const auto &wi : weighted_inits) {
            auto   FObjQ        = wi.init(Q, init_charges);

            // nf_local: number of "light" quarks for this init
            //   ZM (wi.actnf<0): Q-dependent from thresholds
            //   massive (wi.actnf>=0): fixed count of zero-mass quarks
            int    nf_local     = (wi.actnf < 0) ? q2_data[iq].nf : wi.actnf;

            // Sum of EW charges over light quarks
            double sum_ch_light = 0;
            for (int k = 1; k <= nf_local && k - 1 < (int)charges.size(); k++)
                sum_ch_light += charges[k - 1];

            for (int ord = 0; ord < 3; ord++) {
                const auto &C = (ord == 0)   ? FObjQ.C0
                                : (ord == 1) ? FObjQ.C1
                                             : FObjQ.C2;
                if (C.count(1) == 0) continue;

                auto                  &target    = q2_data[iq].channel_ops[ord];

                const auto            &ops_light = C.at(1).GetObjects();
                const apfel::Operator &CNS =
                    ops_light.at(apfel::DISNCBasis::CNS);
                const apfel::Operator &CS = ops_light.at(apfel::DISNCBasis::CS);
                const apfel::Operator &CG_zm =
                    ops_light.at(apfel::DISNCBasis::CG);

                // ── Quark channels ───────────────────────────────────────────
                // C_q = ch_q * CNS + (sum_ch_light/6) * (CS - CNS)
                //
                // For ZM (wi.actnf<0): use actual EW charge for all channels.
                //
                // For massive (wi.actnf>=0): iterate ALL num_quark_channels so
                // the PS contribution (sum_ch_light/6)*(CS-CNS) is summed over
                // all nf_max quark channels and convoluted with the full SIGMA
                // distribution.  Heavy channels (ich >= nf_local) get ch_k=0
                // so they do not contribute the NS term.
                //
                // At NNLO CS-CNS = 6*O22ps, so each of the 5 channels adds
                // (sum_ch_light/6)*6*O22ps ⊗ 2*f = sum_ch_light*O22ps ⊗ 2*f,
                // giving 5*sum_ch_light*O22ps ⊗ 2*f = sum_ch_light*O22ps ⊗
                // 10*f, which matches APFEL++'s (sum_ch_light/6)*CS ⊗ SIGMA
                // (SIGMA=10*f).
                for (int ich = 0; ich < num_quark_channels; ich++) {
                    double ch_k;
                    if (wi.actnf < 0) {
                        // ZM: use actual EW charge for all channels
                        ch_k = (ich < (int)charges.size()) ? charges[ich] : 0.0;
                    } else {
                        // Massive: heavy channels (ich >= nf_local) get no NS
                        // charge
                        ch_k = (ich < nf_local && ich < (int)charges.size())
                                   ? charges[ich]
                                   : 0.0;
                    }
                    auto C_q = ch_k * CNS + (sum_ch_light / 6.0) * (CS - CNS);

                    // Heavy-quark singlet correction: CS non-zero only at NNLO.
                    if (wi.actnf >= 0 && ord == 2) {
                        for (int kh = nf_local + 1; kh <= 6; kh++) {
                            if (C.count(kh) == 0) continue;
                            double M_kh = theory.heavy_quark_masses[kh - 1];
                            if (M_kh < 1e-8) continue;
                            double xi_kh = Q * Q / (M_kh * M_kh);
                            if (xi_kh < theory.mass_ximin ||
                                xi_kh > theory.mass_ximax)
                                continue;
                            const auto &ops_h = C.at(kh).GetObjects();
                            auto it_cs = ops_h.find(apfel::DISNCBasis::CS);
                            if (it_cs == ops_h.end()) continue;
                            double ch_kh  = (kh - 1 < (int)charges.size())
                                                ? charges[kh - 1]
                                                : 0.0;
                            C_q          += (ch_kh / 6.0) * it_cs->second;
                        }
                    }

                    add_to_channel(target, ich, C_q, wi.sign);
                }

                // ── Gluon channel ────────────────────────────────────────────
                // Light contribution: sum_ch_light * CG_zm
                // Heavy contribution (massive only): sum_{k>actnf} ch_k *
                // CG_heavy_k Heavy CG is zero at LO → only access for ord>=1
                // (NLO, NNLO).
                if (has_gluon_channel) {
                    auto CG_total = sum_ch_light * CG_zm;

                    if (wi.actnf >= 0 && ord >= 1) {
                        for (int kh = nf_local + 1; kh <= 6; kh++) {
                            if (C.count(kh) == 0) continue;
                            double M_kh = theory.heavy_quark_masses[kh - 1];
                            if (M_kh < 1e-8) continue;
                            double xi_kh = Q * Q / (M_kh * M_kh);
                            if (xi_kh < theory.mass_ximin ||
                                xi_kh > theory.mass_ximax)
                                continue;
                            const auto &ops_h = C.at(kh).GetObjects();
                            auto it_cg = ops_h.find(apfel::DISNCBasis::CG);
                            if (it_cg == ops_h.end()) continue;
                            double ch_kh  = (kh - 1 < (int)charges.size())
                                                ? charges[kh - 1]
                                                : 0.0;
                            CG_total     += ch_kh * it_cg->second;
                        }
                    }

                    add_to_channel(target, gluon_ich, CG_total, wi.sign);
                }
            }
        }
    }

    // 6. Fill subgrids for each (bin, order, channel) using pre-built operators
    for (std::size_t iord = 0; iord < grid_def.orders.size(); iord++) {
        int alpha_s = grid_def.orders[iord].alpha_s;
        if (alpha_s > 2) {
            std::cerr << "  Warning: skipping order alpha_s^" << alpha_s
                      << " (beyond NNLO)" << std::endl;
            continue;
        }

        for (std::size_t ich = 0; ich < grid_def.channels.size(); ich++) {
            for (std::size_t ibin = 0; ibin < grid_def.bins.size(); ibin++) {
                double              x_lo     = grid_def.bins[ibin].lower.back();
                double              x_hi     = grid_def.bins[ibin].upper.back();
                double              x_center = std::sqrt(x_lo * x_hi);

                std::vector<double> subgrid(nq * nx, 0.0);

                for (std::size_t iq = 0; iq < nq; iq++) {
                    auto it_ord = q2_data[iq].channel_ops.find(alpha_s);
                    if (it_ord == q2_data[iq].channel_ops.end()) continue;

                    auto it_ch = it_ord->second.find((int)ich);
                    if (it_ch == it_ord->second.end()) continue;

                    apfel::Distribution dist = it_ch->second.Evaluate(x_center);
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
