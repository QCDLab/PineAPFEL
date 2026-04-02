#include <apfel/apfelxx.h>
#include <apfel/betaqcd.h>
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
#include <tuple>
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
    int                        n_intermediate       = 3,
    bool                       include_thresholds   = true,
    bool                       use_bin_centers_only = false) {
    std::set<double> q2_set;

    for (const auto &bin : bins) {
        double q2_lo = bin.lower[0];
        double q2_hi = bin.upper[0];
        if (use_bin_centers_only) {
            q2_set.insert(std::sqrt(q2_lo * q2_hi));
            continue;
        }

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

    if (include_thresholds && !q2_set.empty()) {
        double q2_min = *q2_set.begin();
        double q2_max = *q2_set.rbegin();
        for (double thr : thresholds) {
            double q2_thr = thr * thr;
            if (q2_thr > q2_min && q2_thr < q2_max) q2_set.insert(q2_thr);
        }
    }

    return std::vector<double>(q2_set.begin(), q2_set.end());
}

// Evaluate a DoubleOperator at continuous output point (x_c, z_c) and return
// the 2D kernel column w[jx * nj2 + jz] for all joint-grid input nodes (jx,
// jz).
//
// Semantics match APFEL++ `DoubleOperator::operator*=(DoubleDistribution)` for
// a Kronecker input on the joint grid, followed by
// `DoubleDistribution::Evaluate` at (x_c, z_c).  PineAPFEL indexes `jx`, `jz`
// with `GetJointGrid().GetGrid()` sizes (same as `nx_full` / `nz_full` in
// `build_grid_sidis`).
//
// Implementation: for column (jx0, jz0) the joint output before interpolation
// is
//   j(β,δ) = Σ_{α,γ} K(β,δ;α,γ) dj(jsmap1(α), jsmap2(γ))
// with dj = δ_{jx0,jz0}.  Then Evaluate sums L_β(x_c) L_δ(z_c) j(β,δ) over the
// Lagrange support.  Swapping sums yields one pass over (β,δ) in the support
// and (α,γ) in the operator range — O(|support|² · n_α · n_γ) for the full
// column, equivalent to native APFEL++ contraction without building a
// `DoubleDistribution` per joint node.
static std::vector<double> eval_double_op_column(
    const apfel::DoubleOperator       &dop,
    double                             x_c,
    double                             z_c,
    const apfel::LagrangeInterpolator &li1,
    const apfel::LagrangeInterpolator &li2) {
    const apfel::Grid    &g1  = dop.GetFirstGrid();
    const apfel::Grid    &g2  = dop.GetSecondGrid();
    const auto            K   = dop.GetDoubleOperator();

    const std::size_t     nj1 = g1.GetJointGrid().GetGrid().size();
    const std::size_t     nj2 = g2.GetJointGrid().GetGrid().size();
    std::vector<double>   w(nj1 * nj2, 0.0);

    const int             nx1    = g1.GetJointGrid().nx();
    const int             nx2    = g2.GetJointGrid().nx();

    const auto           &jsmap1 = g1.JointToSubMap();
    const auto           &jsmap2 = g2.JointToSubMap();
    const auto           &sjmap1 = g1.SubToJointMap();
    const auto           &sjmap2 = g2.SubToJointMap();

    const apfel::SubGrid &jg1    = g1.GetJointGrid();
    const apfel::SubGrid &jg2    = g2.GetJointGrid();
    const auto            bx     = li1.SumBounds(x_c, jg1);
    const auto            bz     = li2.SumBounds(z_c, jg2);

    for (int beta = bx[0]; beta < bx[1]; beta++) {
        if (beta < 0 || beta >= nx1 || beta >= static_cast<int>(sjmap1.size()))
            continue;
        const double wx = li1.Interpolant(beta, x_c, jg1);
        const auto   m1 = sjmap1[beta];

        for (int delta = bz[0]; delta < bz[1]; delta++) {
            if (delta < 0 || delta >= nx2 ||
                delta >= static_cast<int>(sjmap2.size()))
                continue;
            const double wz   = li2.Interpolant(delta, z_c, jg2);
            const double pref = wx * wz;
            const auto   m2   = sjmap2[delta];

            const auto  &Ksub = K[m1.first][m2.first];
            for (int alpha = m1.second; alpha < g1.GetSubGrid(m1.first).nx();
                 alpha++) {
                const int jx = jsmap1[m1.first][alpha];
                if (jx < 0 || jx >= static_cast<int>(nj1)) continue;
                const int ka = alpha - m1.second;
                if (ka < 0 || ka >= static_cast<int>(Ksub.size(1))) continue;
                const auto &Krow = Ksub(0, ka);
                for (int gamma = m2.second;
                     gamma < g2.GetSubGrid(m2.first).nx();
                     gamma++) {
                    const int jz = jsmap2[m2.first][gamma];
                    if (jz < 0 || jz >= static_cast<int>(nj2)) continue;
                    const int kg = gamma - m2.second;
                    if (kg < 0 || kg >= static_cast<int>(Krow.size(1)))
                        continue;
                    w[static_cast<std::size_t>(jx) * nj2 +
                        static_cast<std::size_t>(jz)] += pref * Krow(0, kg);
                }
            }
        }
    }

    return w;
}

// 8-point Gauss–Legendre nodes/weights on [-1, 1] (∫_{-1}^{1} f ≈ Σ w f(x)).
// Mapped to [z_lo, z_hi] for ∫ G₁(x, z) dz as in APFEL++ BuildSidisG1:
//   Evaluate1(x).Integrate(z_lo, z_hi)
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

// z_out ∈ [z_lo, z_hi]; w_out is the quadrature weight (∫_{z_lo}^{z_hi} f(z) dz
// ≈ Σ w_out f(z_out))
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

// Signed EM charges for quarks 1..6 in APFEL convention (d, u, s, c, b, t)
static constexpr double QCh[6] =
    {-1. / 3., 2. / 3., -1. / 3., 2. / 3., -1. / 3., 2. / 3.};

// Internal implementation: accepts an optional pre-built SidisNNLOObjects.
// When prebuilt_sobj != nullptr, the initialization step is skipped and the
// caller's object is used directly — useful in tests that already own one.
static pineappl_grid *build_grid_sidis(const GridDef &grid_def_in,
    const TheoryCard                                 &theory,
    const OperatorCard                               &op_card,
    const apfel::SidisNNLOObjects                    *prebuilt_sobj = nullptr) {
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

    // 1. Build APFEL++ x/z grids.
    std::vector<apfel::SubGrid> subgrids_x;
    for (const auto &sg : op_card.xgrid)
        subgrids_x.emplace_back(sg.n_knots, sg.x_min, sg.poly_degree);
    std::vector<apfel::SubGrid> subgrids_z = subgrids_x;
    if (!op_card.zgrid.empty()) {
        subgrids_z.clear();
        for (const auto &sg : op_card.zgrid)
            subgrids_z.emplace_back(sg.n_knots, sg.x_min, sg.poly_degree);
    }
    const apfel::Grid gx{subgrids_x};
    const apfel::Grid gz{subgrids_z};

    // 2. Initialize SIDIS coefficient functions (exact NNLO).
    // All orders (LO, NLO, NNLO) use DoubleOperator objects from
    // SidisNNLOObjects.  If the caller already has one, reuse it to
    // avoid a second (expensive) InitializeSidisObjects call.
    using CoeffNNLOPtr = const apfel::DoubleOperator *;

    std::shared_ptr<apfel::SidisNNLOObjects> sobj_owned;
    const apfel::SidisNNLOObjects           *sobj_ptr = prebuilt_sobj;
    if (sobj_ptr == nullptr) {
        sobj_owned =
            std::make_shared<apfel::SidisNNLOObjects>(init_sidis_nnlo(gx,
                gz,
                theory.quark_thresholds,
                op_card.sidis_int_eps));
        sobj_ptr = sobj_owned.get();
    }

    std::function<CoeffNNLOPtr(int, int, Observable, int)> get_coeff =
        [sobj_ptr, polarized](int alpha_s,
            int                   type_id,
            Observable            obs,
            int                   nf) -> CoeffNNLOPtr {
        return select_sidis_coeff_nnlo(*sobj_ptr,
            alpha_s,
            type_id,
            obs,
            nf,
            polarized);
    };

    // 3. Create empty PineAPPL grid
    pineappl_grid           *grid             = create_grid(grid_def);

    // 4. Determine grid nodes
    //
    // APFL++ Lagrange degree-3 grids append a small number of "extension"
    // nodes beyond x = 1 at the upper boundary of the last SubGrid.  These
    // nodes are needed internally by APFL++ for the polynomial basis but
    // have no physical meaning for PDFs/FFs (support is x ≤ 1).  Storing
    // them in the PineAPPL subgrid axis would cause any real PDF library to
    // be called with x > 1 during convolution, which is undefined / an error.
    // Filter them out: keep only nodes with x ≤ 1 and remap the kernel
    // w-array indices accordingly.
    const auto              &x_joint_grid_vec = gx.GetJointGrid().GetGrid();
    const auto              &z_joint_grid_vec = gz.GetJointGrid().GetGrid();
    const std::size_t        nx_full          = x_joint_grid_vec.size();
    const std::size_t        nz_full          = z_joint_grid_vec.size();
    std::vector<double>      x_nodes;
    std::vector<double>      z_nodes;
    std::vector<std::size_t> phys_x_indices;
    std::vector<std::size_t> phys_z_indices;
    x_nodes.reserve(nx_full);
    z_nodes.reserve(nz_full);
    phys_x_indices.reserve(nx_full);
    phys_z_indices.reserve(nz_full);
    for (std::size_t i = 0; i < nx_full; i++) {
        if (x_joint_grid_vec[i] <= 1.0) {
            x_nodes.push_back(x_joint_grid_vec[i]);
            phys_x_indices.push_back(i);
        }
    }
    for (std::size_t i = 0; i < nz_full; i++) {
        if (z_joint_grid_vec[i] <= 1.0) {
            z_nodes.push_back(z_joint_grid_vec[i]);
            phys_z_indices.push_back(i);
        }
    }
    const std::size_t   nx       = x_nodes.size(); // physical x nodes
    const std::size_t   nz       = z_nodes.size(); // physical z nodes

    std::vector<double> q2_nodes = derive_q2_nodes(grid_def.bins,
        theory.quark_thresholds,
        op_card.sidis_q2_n_intermediate,
        op_card.sidis_q2_include_thresholds,
        op_card.sidis_q2_use_bin_centers_only);
    const std::size_t   nq       = q2_nodes.size();

    std::cout << "  Grid nodes: " << nq << " Q^2 x " << nx << " x x " << nz
              << " z points" << std::endl;

    // node_values: [q2_0..q2_{nq-1}, x_0..x_{nx-1}, z_0..z_{nz-1}]
    std::vector<double> node_values;
    node_values.reserve(nq + nx + nz);
    node_values.insert(node_values.end(), q2_nodes.begin(), q2_nodes.end());
    node_values.insert(node_values.end(), x_nodes.begin(), x_nodes.end());
    node_values.insert(node_values.end(), z_nodes.begin(), z_nodes.end());

    std::vector<std::size_t> shape = {nq, nx, nz};

    // 5. Precompute nf per Q² node.
    // IMPORTANT:
    // APFEL++ sidisbuilder uses fixed electromagnetic charge factors QCh2
    // in BuildChannels (not ElectroWeakCharges(Q,false)).
    // To match BuildSidisG1/BuildSidisUnpFT conventions exactly, SIDIS fill
    // must use Q-independent QCh2 weights.
    struct Q2DataSidis {
        int nf;
    };
    std::vector<Q2DataSidis> q2_data(nq);
    for (std::size_t iq = 0; iq < nq; iq++) {
        double Q       = std::sqrt(q2_nodes[iq]);
        q2_data[iq].nf = apfel::NF(Q, theory.quark_thresholds);
    }

    // Lagrange interpolators for NNLO DoubleOperator evaluation.
    const apfel::LagrangeInterpolator li1{gx};
    const apfel::LagrangeInterpolator li2{gz};
    const bool use_bsf_exact = (op_card.sidis_mode == "bsf_exact");

    // 6. Fill subgrids for each (order, channel, bin)
    //
    // Convention: PineAPPL calls the PDF/FF callbacks in LHAPDF convention
    // (returning xf(x)) and divides by x internally before multiplying by the
    // stored subgrid weight.  APFEL++ BuildSidisG1 takes xf inputs and applies
    // an explicit 1/(x·z) factor, so its result is:
    //
    //   G1_BSF(x_c,z_c) = (1/(x_c·z_c)) · Σ_{mn} K_{mn}(x_c,z_c) · xf(x_m) ·
    //   zD(z_n)
    //                    = (1/(x_c·z_c)) · Σ K · x_m·f(x_m) · z_n·D(z_n)
    //
    // PineAPPL with the stored weight w_{mn} computes:
    //   result = Σ αs^n · w_{mn} · f(x_m) · D(z_n)
    //
    // For result == ∫_{z_lo}^{z_hi} G1_BSF(x_c, z) dz (same as
    // BuildSidisG1(...)(Q).Evaluate1(x_c).Integrate(z_lo, z_hi) in
    // sidisbuilder.h and bench_sidis_zurich/sidis_check_zurich.cc) we
    // quadrature-integrate the external z variable:
    //   w_{mn} = Σ_q w_q · K_{mn}(x_c, z_q) · (x_m · z_n) / (x_c · z_q ·
    //   (4π)^n)
    //
    // The (4π)^n factor matches PineAPPL's αs^n against BSF's (αs/4π)^n.
    // The (x_m·z_n)/(x_c·z_ext) factor compensates for the 1/(x_c·z_ext) in
    // BSF.
    for (std::size_t iord = 0; iord < grid_def.orders.size(); iord++) {
        int alpha_s = grid_def.orders[iord].alpha_s;
        if (alpha_s > 2) {
            std::cerr << "  Warning: skipping order alpha_s^" << alpha_s
                      << " (beyond NNLO)" << std::endl;
            continue;
        }
        // Precompute (4π)^alpha_s once per perturbative order.
        const double fourpi_pow =
            std::pow(4.0 * M_PI, static_cast<double>(alpha_s));

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
                double       x_lo         = grid_def.bins[ibin].lower[1];
                double       x_hi         = grid_def.bins[ibin].upper[1];
                double       x_center     = std::sqrt(x_lo * x_hi);

                double       z_lo         = grid_def.bins[ibin].lower[2];
                double       z_hi         = grid_def.bins[ibin].upper[2];
                const bool   z_point_mode = (std::abs(z_hi - z_lo) < 1e-15);

                // When `sidis_q2_use_bin_centers_only` is on, `derive_q2_nodes`
                // inserts exactly one Q² per bin (sqrt(lower·upper)), merged
                // into a global list.  PineAPPL sums the full Q² axis at
                // convolution time, so every slice must be zero except the one
                // matching this bin's center — otherwise each bin picks up
                // every other bin's scale (fatal for multi-bin point grids).
                // With the default tabulation (lo/hi + intermediates), centers
                // need not appear in `q2_nodes`; keep the legacy fill-all-iq
                // path there for backward compatibility.
                const double q2_bin_c = std::sqrt(grid_def.bins[ibin].lower[0] *
                                                  grid_def.bins[ibin].upper[0]);
                const double q2_match_eps =
                    1e-9 * std::max(std::abs(q2_bin_c), 1.0);

                std::vector<double> subgrid(nq * nx * nz, 0.0);
                int                 z_int_lo = 0, z_int_hi = 0;
                std::vector<double> z_int_w;
                if (use_bsf_exact) {
                    if (!z_point_mode) {
                        const apfel::SubGrid &z_joint_sg = gz.GetJointGrid();
                        const auto boundsa2 = li2.SumBounds(z_lo, z_joint_sg);
                        const auto boundsb2 = li2.SumBounds(z_hi, z_joint_sg);
                        z_int_lo            = boundsa2[0];
                        z_int_hi            = boundsb2[1];
                        z_int_w.resize(z_int_hi - z_int_lo);
                        for (int beta = z_int_lo; beta < z_int_hi; beta++) {
                            z_int_w[beta - z_int_lo] = li2.IntInterpolant(beta,
                                z_lo,
                                z_hi,
                                z_joint_sg);
                        }
                    }
                }

                // NOTE: Skip bins where `x_center` is outside APFEL++
                // interpolation range. External z is handled by quadrature on
                // [z_lo, z_hi]; individual points must lie on the joint grid.
                if (x_center >= x_nodes.front() && x_center <= x_nodes.back()) {
                    for (std::size_t iq = 0; iq < nq; iq++) {
                        if (op_card.sidis_q2_use_bin_centers_only &&
                            std::abs(q2_nodes[iq] - q2_bin_c) > q2_match_eps)
                            continue;
                        int                     nf          = q2_data[iq].nf;
                        static constexpr double QCh2[6]     = {1. / 9.,
                                4. / 9.,
                                1. / 9.,
                                4. / 9.,
                                1. / 9.,
                                4. / 9.};

                        // Compute fill weight; skip if required quarks
                        // inactive. For gg/ps the weight is sumch2 = Σ_q e_q²;
                        // for qpq3 it is QCh[a-1]*QCh[b-1] (may be negative).
                        double                  fill_weight = 0.0;
                        switch (info.type) {
                        case SidisChannelInfo::Type::nn:
                        case SidisChannelInfo::Type::gq:
                        case SidisChannelInfo::Type::qg:
                            if (info.quark_a > nf) continue;
                            fill_weight = QCh2[info.quark_a - 1];
                            break;
                        case SidisChannelInfo::Type::gg:
                            for (int q = 0; q < nf; q++) fill_weight += QCh2[q];
                            break;
                        case SidisChannelInfo::Type::ps:
                            if (info.quark_a > nf) continue;
                            for (int q = 0; q < nf; q++) fill_weight += QCh2[q];
                            break;
                        case SidisChannelInfo::Type::qbq:
                        case SidisChannelInfo::Type::qpq1:
                        case SidisChannelInfo::Type::qpq2:
                            if (info.quark_a > nf) continue;
                            fill_weight = QCh2[info.quark_a - 1];
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

                        if (use_bsf_exact) {
                            if (z_point_mode) {
                                const double z_out = z_lo;
                                if (z_out < z_nodes.front() ||
                                    z_out > z_nodes.back())
                                    continue;
                                auto         w = eval_double_op_column(*coeff,
                                    x_center,
                                    z_out,
                                    li1,
                                    li2);
                                const double xz_denom =
                                    x_center * z_out * fourpi_pow;
                                for (std::size_t iix = 0; iix < nx; iix++) {
                                    const std::size_t ix = phys_x_indices[iix];
                                    const double      x_node = x_nodes[iix];
                                    for (std::size_t iiz = 0; iiz < nz; iiz++) {
                                        const std::size_t iz =
                                            phys_z_indices[iiz];
                                        const double z_node = z_nodes[iiz];
                                        subgrid[iq * nx * nz + iix * nz +
                                                iiz] +=
                                            fill_weight * w[ix * nz_full + iz] *
                                            (x_node * z_node) / xz_denom;
                                    }
                                }
                            } else {
                                const auto &z_joint_nodes =
                                    gz.GetJointGrid().GetGrid();
                                for (int beta = z_int_lo; beta < z_int_hi;
                                     beta++) {
                                    const double z_out = z_joint_nodes[beta];
                                    if (z_out < z_nodes.front() ||
                                        z_out > z_nodes.back())
                                        continue;

                                    auto w = eval_double_op_column(*coeff,
                                        x_center,
                                        z_out,
                                        li1,
                                        li2);

                                    // APFEL++ does z-integration as
                                    // IntInterpolant weights times nodal
                                    // values.
                                    const double xz_denom =
                                        x_center * z_out * fourpi_pow;
                                    const double int_w =
                                        z_int_w[beta - z_int_lo];
                                    for (std::size_t iix = 0; iix < nx; iix++) {
                                        const std::size_t ix =
                                            phys_x_indices[iix];
                                        const double x_node = x_nodes[iix];
                                        for (std::size_t iiz = 0; iiz < nz;
                                             iiz++) {
                                            const std::size_t iz =
                                                phys_z_indices[iiz];
                                            const double z_node = z_nodes[iiz];
                                            subgrid[iq * nx * nz + iix * nz +
                                                    iiz] +=
                                                fill_weight * int_w *
                                                w[ix * nz_full + iz] *
                                                (x_node * z_node) / xz_denom;
                                        }
                                    }
                                }
                            }
                        } else {
                            if (z_point_mode) {
                                const double z_ext = z_lo;
                                if (z_ext < z_nodes.front() ||
                                    z_ext > z_nodes.back())
                                    continue;
                                auto         w = eval_double_op_column(*coeff,
                                    x_center,
                                    z_ext,
                                    li1,
                                    li2);
                                const double xz_denom =
                                    x_center * z_ext * fourpi_pow;
                                for (std::size_t iix = 0; iix < nx; iix++) {
                                    const std::size_t ix = phys_x_indices[iix];
                                    const double      x_node = x_nodes[iix];
                                    for (std::size_t iiz = 0; iiz < nz; iiz++) {
                                        const std::size_t iz =
                                            phys_z_indices[iiz];
                                        const double z_node = z_nodes[iiz];
                                        subgrid[iq * nx * nz + iix * nz +
                                                iiz] +=
                                            fill_weight * w[ix * nz_full + iz] *
                                            (x_node * z_node) / xz_denom;
                                    }
                                }
                            } else {
                                const int z_nsub =
                                    op_card.sidis_z_quad_subdivisions;
                                for (int isub = 0; isub < z_nsub; isub++) {
                                    const double z_a = z_lo + (z_hi - z_lo) *
                                                                  (double)isub /
                                                                  z_nsub;
                                    const double z_b =
                                        z_lo + (z_hi - z_lo) *
                                                   (double)(isub + 1) / z_nsub;
                                    for (int izq = 0; izq < kSidisZQuadN;
                                         izq++) {
                                        double z_ext = 0.0, dz_w = 0.0;
                                        sidis_z_quadrature_point(z_a,
                                            z_b,
                                            izq,
                                            z_ext,
                                            dz_w);
                                        if (z_ext < z_nodes.front() ||
                                            z_ext > z_nodes.back())
                                            continue;

                                        auto w = eval_double_op_column(*coeff,
                                            x_center,
                                            z_ext,
                                            li1,
                                            li2);

                                        const double xz_denom =
                                            x_center * z_ext * fourpi_pow;
                                        for (std::size_t iix = 0; iix < nx;
                                             iix++) {
                                            const std::size_t ix =
                                                phys_x_indices[iix];
                                            const double x_node = x_nodes[iix];
                                            for (std::size_t iiz = 0; iiz < nz;
                                                 iiz++) {
                                                const std::size_t iz =
                                                    phys_z_indices[iiz];
                                                const double z_node =
                                                    z_nodes[iiz];
                                                subgrid[iq * nx * nz +
                                                        iix * nz + iiz] +=
                                                    fill_weight * dz_w *
                                                    w[ix * nz_full + iz] *
                                                    (x_node * z_node) /
                                                    xz_denom;
                                            }
                                        }
                                    }
                                }
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

    // 4. Determine grid nodes — filter out APFL++ Lagrange extension nodes
    //    beyond x = 1 (same rationale as in build_grid_sidis).
    const auto         &joint_grid_vec = g.GetJointGrid().GetGrid();
    const std::size_t   nx_full        = joint_grid_vec.size();
    std::vector<double> x_nodes;
    std::vector<std::size_t> phys_x_indices;
    x_nodes.reserve(nx_full);
    phys_x_indices.reserve(nx_full);
    for (std::size_t i = 0; i < nx_full; i++) {
        if (joint_grid_vec[i] <= 1.0) {
            x_nodes.push_back(joint_grid_vec[i]);
            phys_x_indices.push_back(i);
        }
    }
    const std::size_t   nx       = x_nodes.size();

    // Pointwise design: one Q² node per bin at the geometric bin centre.
    // Collect unique Q²_c values to avoid redundant APFEL++ coefficient
    // function calls when multiple bins share the same Q²_c.
    std::vector<double> bin_q2c(grid_def.bins.size());
    std::set<double>    unique_q2c_set;
    for (std::size_t ibin = 0; ibin < grid_def.bins.size(); ibin++) {
        bin_q2c[ibin] = std::sqrt(grid_def.bins[ibin].lower[0] *
                                  grid_def.bins[ibin].upper[0]);
        unique_q2c_set.insert(bin_q2c[ibin]);
    }
    std::vector<double> unique_q2c(unique_q2c_set.begin(), unique_q2c_set.end());

    std::cout << "  Pointwise grid: " << unique_q2c.size()
              << " unique Q² centres x " << nx << " x/z points" << std::endl;

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
    //    channel_ops[{alpha_s, log_xir, log_xif}][channel_idx] accumulates
    //    the combined operator from all weighted initializers.  Gluon and
    //    quark channels are handled separately to correctly incorporate
    //    heavy-quark contributions in the massive scheme.
    //    Scale-log orders ({alpha_s, >=1, 0} or {alpha_s, 0, >=1}) are
    //    derived from the central-scale operators after the main accumulation.
    using OrdKey = std::tuple<int, int, int>; // (alpha_s, log_xir, log_xif)
    struct Q2Data {
        int                                              nf;
        // 6-entry EW charges (NC) or nf-entry CKM (CC):
        std::vector<double>                              charges;
        std::map<OrdKey, std::map<int, apfel::Operator>> channel_ops;
    };

    std::map<double, Q2Data> q2_data_map;

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

    for (double q2_c : unique_q2c) {
        Q2Data &qd     = q2_data_map[q2_c];
        double  Q      = std::sqrt(q2_c);
        qd.nf          = apfel::NF(Q, theory.quark_thresholds);

        // Compute per-quark charges and initializer charges
        std::vector<double> init_charges;
        if (is_cc) {
            int                 nf = qd.nf;
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
            qd.charges   = weights;
            init_charges = theory.ckm;
        } else {
            // ElectroWeakCharges always returns 6 entries
            qd.charges   = apfel::ElectroWeakCharges(Q, timelike);
            init_charges = qd.charges;
        }

        const std::vector<double> &charges = qd.charges;

        for (const auto &wi : weighted_inits) {
            auto   FObjQ        = wi.init(Q, init_charges);

            // nf_local: number of "light" quarks for this init
            //   ZM (wi.actnf<0): Q-dependent from thresholds
            //   massive (wi.actnf>=0): fixed count of zero-mass quarks
            int    nf_local     = (wi.actnf < 0) ? qd.nf : wi.actnf;

            // Sum of EW charges over light quarks
            double sum_ch_light = 0;
            for (int k = 1; k <= nf_local && k - 1 < (int)charges.size(); k++)
                sum_ch_light += charges[k - 1];

            for (int ord = 0; ord < 3; ord++) {
                const auto &C = (ord == 0)   ? FObjQ.C0
                                : (ord == 1) ? FObjQ.C1
                                             : FObjQ.C2;
                if (C.count(1) == 0) continue;

                // Scale by (1/(4π))^ord so that when PineAPPL multiplies by
                // α_s^ord the result matches APFEL++ BSF which uses
                // (α_s/4π)^ord.
                const double pert_norm = std::pow(1.0 / apfel::FourPi, ord);
                const OrdKey key_ord{ord, 0, 0};
                auto        &target    = qd.channel_ops[key_ord];

                const auto  &ops_light = C.at(1).GetObjects();
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

                    add_to_channel(target, ich, C_q * pert_norm, wi.sign);
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

                    add_to_channel(target,
                        gluon_ich,
                        CG_total * pert_norm,
                        wi.sign);
                }
            }
        }

        // Renormalization-scale logs derived from central-scale operators.
        // Populate channel_ops[{n, m, 0}] for m>0 from the accumulated
        // central-scale channel_ops[{n-m, 0, 0}] entries.
        //
        // Coefficients follow from d/d(tR) [C(ξ_R)] with tR = ln(ξ_R²):
        //   {1,1,0}: b₀ · C0 · (1/4π)^1
        //   {2,1,0}: b₀ · C1 · (1/4π)^2
        //   {2,2,0}: (b₀²/2) · C0 · (1/4π)^2
        // where each C_n is already stored in channel_ops[{n,0,0}] with its
        // (1/4π)^n factor absorbed, so the extra factor per log-slot is just
        // (1/4π)^{alpha_s} / (1/4π)^{base_ord} = (1/4π)^{log_xir}.
        // Only computed when the requested order set contains log_xir > 0.
        {
            bool need_xir = false;
            for (const auto &o : grid_def.orders)
                if (o.log_xir > 0) {
                    need_xir = true;
                    break;
                }

            if (need_xir) {
                const double b0     = apfel::beta0qcd(qd.nf);
                const double inv4pi = 1.0 / apfel::FourPi;
                auto        &cops   = qd.channel_ops;

                auto         it0    = cops.find(OrdKey{0, 0, 0});
                auto         it1    = cops.find(OrdKey{1, 0, 0});

                // {1,1,0}: b₀ * (1/4π) * C0_channels
                if (it0 != cops.end()) {
                    const double fac = b0 * inv4pi;
                    for (const auto &[ich, op] : it0->second)
                        add_to_channel(cops[OrdKey{1, 1, 0}],
                            ich,
                            op * fac,
                            1.0);
                }
                // {2,1,0}: b₀ * (1/4π) * C1_channels
                if (it1 != cops.end()) {
                    const double fac = b0 * inv4pi;
                    for (const auto &[ich, op] : it1->second)
                        add_to_channel(cops[OrdKey{2, 1, 0}],
                            ich,
                            op * fac,
                            1.0);
                }
                // {2,2,0}: (b₀²/2) * (1/4π)² * C0_channels
                if (it0 != cops.end()) {
                    const double fac = 0.5 * b0 * b0 * inv4pi * inv4pi;
                    for (const auto &[ich, op] : it0->second)
                        add_to_channel(cops[OrdKey{2, 2, 0}],
                            ich,
                            op * fac,
                            1.0);
                }
            }
        }
    }

    // 6. Fill subgrids for each (bin, order, channel) using pre-built operators
    for (std::size_t iord = 0; iord < grid_def.orders.size(); iord++) {
        const auto &odef    = grid_def.orders[iord];
        int         alpha_s = odef.alpha_s;
        if (alpha_s > 2) {
            std::cerr << "  Warning: skipping order alpha_s^" << alpha_s
                      << " (beyond NNLO)" << std::endl;
            continue;
        }
        const OrdKey ord_key{(int)odef.alpha_s,
            (int)odef.log_xir,
            (int)odef.log_xif};

        for (std::size_t ich = 0; ich < grid_def.channels.size(); ich++) {
            for (std::size_t ibin = 0; ibin < grid_def.bins.size(); ibin++) {
                double x_lo     = grid_def.bins[ibin].lower.back();
                double x_hi     = grid_def.bins[ibin].upper.back();
                double x_center = std::sqrt(x_lo * x_hi);

                // Per-bin 1-node Q² subgrid at the geometric bin centre.
                double              q2_c = bin_q2c[ibin];
                std::vector<double> node_vals_bin;
                node_vals_bin.reserve(1 + nx);
                node_vals_bin.push_back(q2_c);
                node_vals_bin.insert(node_vals_bin.end(),
                    x_nodes.begin(), x_nodes.end());
                std::vector<std::size_t> shape_bin = {1, nx};

                std::vector<double> subgrid(nx, 0.0);

                // NOTE: Skip bins outside the APFEL++ interpolation range; the
                // operator Evaluate() has undefined behaviour for x outside
                // [x_nodes.front(), x_nodes.back()].
                if (x_center >= x_nodes.front() && x_center <= x_nodes.back()) {
                    auto it_q2 = q2_data_map.find(q2_c);
                    if (it_q2 != q2_data_map.end()) {
                        auto it_ord = it_q2->second.channel_ops.find(ord_key);
                        if (it_ord != it_q2->second.channel_ops.end()) {
                            auto it_ch = it_ord->second.find((int)ich);
                            if (it_ch != it_ord->second.end()) {
                                apfel::Distribution dist =
                                    it_ch->second.Evaluate(x_center);
                                const std::vector<double> &vals =
                                    dist.GetDistributionJointGrid();

                                // NOTE: APFEL++: K_j satisfy (C⊗f)(x)=Σ K_j*f(x_j)
                                // PineAPPL divides xfx(x_j) by x_j, so store K_j*x_j.
                                // vals is indexed over the full joint grid (nx_full).
                                for (std::size_t iix = 0; iix < nx; iix++) {
                                    const std::size_t ix = phys_x_indices[iix];
                                    if (ix < vals.size())
                                        subgrid[iix] = vals[ix] * x_nodes[iix];
                                }
                            }
                        }
                    }
                }

                set_subgrid(grid, ibin, iord, ich, node_vals_bin, subgrid, shape_bin);
            }
        }
    }

    std::cout << "Grid filled successfully." << std::endl;
    return grid;
}

pineappl_grid *build_grid(const GridDef &grid_def,
    const TheoryCard                    &theory,
    const OperatorCard                  &op_card,
    const apfel::SidisNNLOObjects       &prebuilt_sobj) {
    if (grid_def.process == ProcessType::SIDIS)
        return build_grid_sidis(grid_def, theory, op_card, &prebuilt_sobj);
    // Non-SIDIS: prebuilt_sobj is irrelevant; delegate to the standard
    // overload.
    return build_grid(grid_def, theory, op_card);
}

} // namespace pineapfel
