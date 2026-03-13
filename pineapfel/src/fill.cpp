#include <apfel/apfelxx.h>
#include <pineapfel/fill.h>
#include <pineapfel/sidis_api.h>

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

// Evaluate a DoubleOperator at continuous output point (x_c, z_c) and return
// the 2D kernel column w[ix * nz + iz] for all joint-grid input nodes (ix, iz).
// This is the 2D analogue of Operator::Evaluate(x_c) → Distribution for
// PineAPPL subgrid filling: w[ix][iz] = K(x_c, z_c; x_{ix}, z_{iz}).
//
// The shift-invariant kernel K_raw[ig1][ig2](0, beta_x)(0, beta_z) stores the
// coupling for x-separation beta_x and z-separation beta_z within subgrid pair
// (ig1, ig2).  The output-point Lagrange weights are used to combine
// contributions from grid nodes near (x_c, z_c).
static std::vector<double> eval_double_op_column(
    const apfel::DoubleOperator       &dop,
    double                             x_c,
    double                             z_c,
    const apfel::LagrangeInterpolator &li1,
    const apfel::LagrangeInterpolator &li2) {
    const apfel::Grid &g1  = dop.GetFirstGrid();
    const apfel::Grid &g2  = dop.GetSecondGrid();
    const auto        &K   = dop.GetDoubleOperator();

    const int          ng1 = g1.nGrids();
    const int          ng2 = g2.nGrids();

    // Total joint-grid sizes
    int                nx  = 0;
    for (int ig = 0; ig < ng1; ig++) nx += g1.GetSubGrid(ig).nx();
    int nz = 0;
    for (int ig = 0; ig < ng2; ig++) nz += g2.GetSubGrid(ig).nx();

    std::vector<double> w(nx * nz, 0.0);

    // Find subgrid ig1 that is the densest one still containing x_c
    int                 ig1_xc = ng1 - 1;
    for (int ig1 = 0; ig1 < ng1; ig1++)
        if (g1.GetSubGrid(ig1).xMin() > x_c) {
            ig1_xc = ig1 - 1;
            break;
        }

    // Find subgrid ig2 containing z_c
    int ig2_zc = ng2 - 1;
    for (int ig2 = 0; ig2 < ng2; ig2++)
        if (g2.GetSubGrid(ig2).xMin() > z_c) {
            ig2_zc = ig2 - 1;
            break;
        }

    const apfel::SubGrid &sg1    = g1.GetSubGrid(ig1_xc);
    const apfel::SubGrid &sg2    = g2.GetSubGrid(ig2_zc);
    const int             nx_sub = sg1.nx();
    const int             nz_sub = sg2.nx();

    // Lagrange weights for x_c within sg1
    const auto            bx     = li1.SumBounds(x_c, sg1);
    const int             lo_x = bx[0], hi_x = bx[1];
    std::vector<double>   W_x(hi_x - lo_x);
    for (int b = 0; b < hi_x - lo_x; b++)
        W_x[b] = li1.Interpolant(lo_x + b, x_c, sg1);

    // Lagrange weights for z_c within sg2
    const auto          bz   = li2.SumBounds(z_c, sg2);
    const int           lo_z = bz[0], hi_z = bz[1];
    std::vector<double> W_z(hi_z - lo_z);
    for (int b = 0; b < hi_z - lo_z; b++)
        W_z[b] = li2.Interpolant(lo_z + b, z_c, sg2);

    // Joint-grid offset for ig1_xc and ig2_zc
    int x_jg_off = 0;
    for (int ig = 0; ig < ig1_xc; ig++) x_jg_off += g1.GetSubGrid(ig).nx();
    int z_jg_off = 0;
    for (int ig = 0; ig < ig2_zc; ig++) z_jg_off += g2.GetSubGrid(ig).nx();

    const auto &K_sub = K[ig1_xc][ig2_zc];

    // For each output node (out_x, out_z) near (x_c, z_c), accumulate
    // kernel weights for all connected input nodes (out_x+alpha, out_z+gamma).
    for (int b_ox = 0; b_ox < hi_x - lo_x; b_ox++) {
        const int    out_x = lo_x + b_ox;
        const double wx    = W_x[b_ox];

        for (int b_oz = 0; b_oz < hi_z - lo_z; b_oz++) {
            const int    out_z = lo_z + b_oz;
            const double wxz   = wx * W_z[b_oz];

            for (int alpha = 0; alpha <= nx_sub - out_x - 1; alpha++) {
                const int inp_x = x_jg_off + out_x + alpha;
                if (inp_x >= nx) break;

                for (int gamma = 0; gamma <= nz_sub - out_z - 1; gamma++) {
                    const int inp_z = z_jg_off + out_z + gamma;
                    if (inp_z >= nz) break;

                    w[inp_x * nz + inp_z] += wxz * K_sub(0, alpha)(0, gamma);
                }
            }
        }
    }

    return w;
}

// Select a DoubleOperator from SidisNNLOObjects for a given
// (alpha_s, type_id, observable, nf, polarized).
// type_id: 0=nn/ns, 1=gq, 2=qg, 3=gg, 4=ps, 5=qbq, 6=qpq1, 7=qpq2, 8=qpq3
// Returns nullptr when no operator exists (e.g., FL at LO, pure-NNLO channels
// at NLO).
static const apfel::DoubleOperator *select_sidis_coeff_nnlo(
    const apfel::SidisNNLOObjects &obj,
    int                            alpha_s,
    int                            type_id,
    Observable                     observable,
    int                            nf,
    bool                           polarized) {
    const std::string snf = "_nf" + std::to_string(nf);
    const auto       &map =
        polarized ? obj.G1 : (observable == Observable::FL ? obj.FL : obj.FT);

    auto find = [&](const std::string &key) -> const apfel::DoubleOperator * {
        auto it = map.find(key);
        return it != map.end() ? &it->second : nullptr;
    };

    if (polarized) {
        // G1 polarized
        if (alpha_s == 0) {
            return type_id == 0 ? find("DoubleIdentity") : nullptr;
        } else if (alpha_s == 1) {
            switch (type_id) {
            case 0 : return find("DC1Q2Q");
            case 1 : return find("DC1Q2G");
            case 2 : return find("DC1G2Q");
            default: return nullptr;
            }
        } else {
            switch (type_id) {
            case 0 : return find("DC2Q2QNS" + snf);
            case 1 : return find("DC2Q2G" + snf);
            case 2 : return find("DC2G2Q" + snf);
            case 3 : return find("DC2G2G");
            case 4 : return find("DC2Q2QPS");
            case 5 : return find("DC2Q2QB");
            case 6 : return find("DC2Q2QP1");
            case 7 : return find("DC2Q2QP2");
            case 8 : return find("DC2Q2QP3");
            default: return nullptr;
            }
        }
    } else if (observable == Observable::F2) {
        // FT unpolarized
        if (alpha_s == 0) {
            return type_id == 0 ? find("DoubleIdentity") : nullptr;
        } else if (alpha_s == 1) {
            switch (type_id) {
            case 0 : return find("C1TQ2Q");
            case 1 : return find("C1TQ2G");
            case 2 : return find("C1TG2Q");
            default: return nullptr;
            }
        } else {
            switch (type_id) {
            case 0 : return find("C2TQ2QNS" + snf);
            case 1 : return find("C2TQ2G");
            case 2 : return find("C2TG2Q");
            case 3 : return find("C2TG2G");
            case 4 : return find("C2TQ2QPS");
            case 5 : return find("C2TQ2QB");
            case 6 : return find("C2TQ2QP1");
            case 7 : return find("C2TQ2QP2");
            case 8 : return find("C2TQ2QP3");
            default: return nullptr;
            }
        }
    } else {
        // FL unpolarized — no LO contribution
        if (alpha_s == 0) return nullptr;
        if (alpha_s == 1) {
            switch (type_id) {
            case 0 : return find("C1LQ2Q");
            case 1 : return find("C1LQ2G");
            case 2 : return find("C1LG2Q");
            default: return nullptr;
            }
        } else {
            switch (type_id) {
            case 0 : return find("C2LQ2QNS" + snf);
            case 1 : return find("C2LQ2G");
            case 2 : return find("C2LG2Q");
            case 3 : return find("C2LG2G");
            case 4 : return find("C2LQ2QPS");
            case 5 : return find("C2LQ2QB");
            case 6 : return find("C2LQ2QP1");
            case 7 : return find("C2LQ2QP2");
            case 8 : return find("C2LQ2QP3");
            default: return nullptr;
            }
        }
    }
    return nullptr;
}

// Describes the type and quark index(es) of a SIDIS channel for fill routing.
struct SidisChannelInfo {
    enum class Type { nn, gq, qg, gg, ps, qbq, qpq1, qpq2, qpq3 } type;
    int quark_a = -1; // primary quark (1-based): active quark for
                      // nn/gq/qg/ps/qbq/qpq1/qpq2
    int quark_b = -1; // secondary quark (1-based): target for qpq3
};

// Signed EM charges for quarks 1..6 (u, d, s, c, b, t)
static constexpr double QCh[6] =
    {2. / 3., -1. / 3., -1. / 3., 2. / 3., -1. / 3., 2. / 3.};

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

    // Build SidisChannelInfo in the same order as derive_channels.
    std::vector<SidisChannelInfo> channel_infos;
    channel_infos.reserve(grid_def.channels.size());
    // 3*nf_max existing channels: nn, gq, qg per quark
    for (int q = 1; q <= nf_max; q++) {
        channel_infos.push_back({SidisChannelInfo::Type::nn, q, -1});
        channel_infos.push_back({SidisChannelInfo::Type::gq, q, -1});
        channel_infos.push_back({SidisChannelInfo::Type::qg, q, -1});
    }
    // gg
    channel_infos.push_back({SidisChannelInfo::Type::gg, -1, -1});
    // ps per quark
    for (int q = 1; q <= nf_max; q++)
        channel_infos.push_back({SidisChannelInfo::Type::ps, q, -1});
    // qbq per quark
    for (int q = 1; q <= nf_max; q++)
        channel_infos.push_back({SidisChannelInfo::Type::qbq, q, -1});
    // qpq1 per source j (only when nf_max >= 2, matching derive_channels guard)
    for (int j = 1; j <= nf_max; j++) {
        if (nf_max < 2) continue;
        channel_infos.push_back({SidisChannelInfo::Type::qpq1, j, -1});
    }
    // qpq2 per target i
    for (int i = 1; i <= nf_max; i++) {
        if (nf_max < 2) continue;
        channel_infos.push_back({SidisChannelInfo::Type::qpq2, i, -1});
    }
    // qpq3 per ordered pair (a,b) with a≠b
    for (int a = 1; a <= nf_max; a++)
        for (int b = 1; b <= nf_max; b++)
            if (a != b)
                channel_infos.push_back({SidisChannelInfo::Type::qpq3, a, b});

    if (channel_infos.size() != grid_def.channels.size())
        throw std::runtime_error("build_grid_sidis: channel_infos size "
                                 "mismatch with derived channels");

    // 1. Build APFEL++ x-space grid (same grid for x and z)
    std::vector<apfel::SubGrid> subgrids;
    for (const auto &sg : op_card.xgrid)
        subgrids.emplace_back(sg.n_knots, sg.x_min, sg.poly_degree);
    const apfel::Grid g{subgrids};

    // 2. Initialize SIDIS coefficient functions (exact NNLO).
    // All orders (LO, NLO, NNLO) use DoubleOperator objects from
    // SidisNNLOObjects.
    using CoeffNNLOPtr = const apfel::DoubleOperator *;

    auto sobj          = std::make_shared<apfel::SidisNNLOObjects>(
        init_sidis_nnlo(g, theory.quark_thresholds, op_card.sidis_int_eps));
    std::function<CoeffNNLOPtr(int, int, Observable, int)> get_coeff =
        [sobj, polarized](int alpha_s,
            int               type_id,
            Observable        obs,
            int               nf) -> CoeffNNLOPtr {
        return select_sidis_coeff_nnlo(*sobj,
            alpha_s,
            type_id,
            obs,
            nf,
            polarized);
    };

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

    // Lagrange interpolators for NNLO DoubleOperator evaluation.
    // Both x and z use the same grid g.
    const apfel::LagrangeInterpolator li1{g};
    const apfel::LagrangeInterpolator li2{g};

    // 6. Fill subgrids for each (order, channel, bin)
    for (std::size_t iord = 0; iord < grid_def.orders.size(); iord++) {
        int alpha_s = grid_def.orders[iord].alpha_s;
        if (alpha_s > 2) {
            std::cerr << "  Warning: skipping order alpha_s^" << alpha_s
                      << " (beyond NNLO)" << std::endl;
            continue;
        }

        for (std::size_t ich = 0; ich < grid_def.channels.size(); ich++) {
            const SidisChannelInfo &info    = channel_infos[ich];

            // Map channel type to unified type_id (0..8).
            int                     type_id = -1;
            switch (info.type) {
            case SidisChannelInfo::Type::nn  : type_id = 0; break;
            case SidisChannelInfo::Type::gq  : type_id = 1; break;
            case SidisChannelInfo::Type::qg  : type_id = 2; break;
            case SidisChannelInfo::Type::gg  : type_id = 3; break;
            case SidisChannelInfo::Type::ps  : type_id = 4; break;
            case SidisChannelInfo::Type::qbq : type_id = 5; break;
            case SidisChannelInfo::Type::qpq1: type_id = 6; break;
            case SidisChannelInfo::Type::qpq2: type_id = 7; break;
            case SidisChannelInfo::Type::qpq3: type_id = 8; break;
            }
            if (type_id < 0) continue;

            for (std::size_t ibin = 0; ibin < grid_def.bins.size(); ibin++) {
                // Bin centers: dim 0=Q², dim 1=x, dim 2=z
                double              x_lo     = grid_def.bins[ibin].lower[1];
                double              x_hi     = grid_def.bins[ibin].upper[1];
                double              x_center = std::sqrt(x_lo * x_hi);

                double              z_lo     = grid_def.bins[ibin].lower[2];
                double              z_hi     = grid_def.bins[ibin].upper[2];
                double              z_center = std::sqrt(z_lo * z_hi);

                std::vector<double> subgrid(nq * nx * nz, 0.0);

                // NOTE: Skip bins where `x_center` or `z_center` are outside
                // APFEL++ interpolation range to avoid undefined behaviour in
                // Operator::Evaluate().  The subgrid stays zero, consistent
                // with PineAPPL returning zero for convolutions below x_min.
                if (x_center >= x_nodes.front() && x_center <= x_nodes.back() &&
                    z_center >= x_nodes.front() && z_center <= x_nodes.back()) {
                    for (std::size_t iq = 0; iq < nq; iq++) {
                        int         nf          = q2_data[iq].nf;
                        const auto &charges     = q2_data[iq].charges;

                        // Compute fill weight; skip if required quarks
                        // inactive. For gg/ps the weight is sumch2 = Σ_q e_q²;
                        // for qpq3 it is QCh[a-1]*QCh[b-1] (may be negative).
                        double      fill_weight = 0.0;
                        switch (info.type) {
                        case SidisChannelInfo::Type::nn:
                        case SidisChannelInfo::Type::gq:
                        case SidisChannelInfo::Type::qg:
                            if (info.quark_a > nf) continue;
                            fill_weight = charges[info.quark_a - 1];
                            break;
                        case SidisChannelInfo::Type::gg:
                            for (int q = 0; q < nf; q++)
                                fill_weight += charges[q];
                            break;
                        case SidisChannelInfo::Type::ps:
                            if (info.quark_a > nf) continue;
                            for (int q = 0; q < nf; q++)
                                fill_weight += charges[q];
                            break;
                        case SidisChannelInfo::Type::qbq:
                        case SidisChannelInfo::Type::qpq1:
                        case SidisChannelInfo::Type::qpq2:
                            if (info.quark_a > nf) continue;
                            fill_weight = charges[info.quark_a - 1];
                            break;
                        case SidisChannelInfo::Type::qpq3:
                            if (info.quark_a > nf || info.quark_b > nf)
                                continue;
                            fill_weight =
                                QCh[info.quark_a - 1] * QCh[info.quark_b - 1];
                            break;
                        }

                        if (fill_weight == 0.0) continue;

                        // All orders use DoubleOperator (exact NNLO API).
                        const auto *coeff = get_coeff(alpha_s,
                            type_id,
                            grid_def.observable,
                            nf);
                        if (coeff == nullptr) continue;

                        auto w = eval_double_op_column(*coeff,
                            x_center,
                            z_center,
                            li1,
                            li2);

                        for (std::size_t ix = 0; ix < nx; ix++) {
                            for (std::size_t iz = 0; iz < nz; iz++) {
                                subgrid[iq * nx * nz + ix * nz + iz] +=
                                    fill_weight * w[ix * nz + iz];
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

                // NOTE: Skip bins outside the APFEL++ interpolation range; the
                // operator Evaluate() has undefined behaviour for x outside
                // [x_nodes.front(), x_nodes.back()].
                if (x_center >= x_nodes.front() && x_center <= x_nodes.back()) {
                    for (std::size_t iq = 0; iq < nq; iq++) {
                        auto it_ord = q2_data[iq].channel_ops.find(alpha_s);
                        if (it_ord == q2_data[iq].channel_ops.end()) continue;

                        auto it_ch = it_ord->second.find((int)ich);
                        if (it_ch == it_ord->second.end()) continue;

                        apfel::Distribution dist =
                            it_ch->second.Evaluate(x_center);
                        const std::vector<double> &vals =
                            dist.GetDistributionJointGrid();

                        for (std::size_t ix = 0; ix < nx && ix < vals.size();
                             ix++)
                            subgrid[iq * nx + ix] = vals[ix];
                    }
                }

                set_subgrid(grid, ibin, iord, ich, node_values, subgrid, shape);
            }
        }
    }

    std::cout << "Grid filled successfully." << std::endl;
    return grid;
}

} // namespace pineapfel
