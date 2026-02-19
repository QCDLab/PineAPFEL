#include <grid_gen.h>
#include <yaml-cpp/yaml.h>

#include <numeric>
#include <stdexcept>

namespace pineapfel {

std::vector<ChannelDef> derive_channels(Observable observable, int nf_max) {
    std::vector<ChannelDef> channels;

    bool is_f3 = (observable == Observable::F3);

    // One channel per active quark flavor: (q + qbar) for F2/FL, (q - qbar) for F3
    for (int q = 1; q <= nf_max; q++) {
        ChannelDef ch;
        ch.pid_combinations = {{q}, {-q}};
        ch.factors = is_f3 ? std::vector<double>{1.0, -1.0}
                           : std::vector<double>{1.0, 1.0};
        channels.push_back(std::move(ch));
    }

    // Gluon channel only for F2/FL (CG=0 for F3 at all orders)
    if (!is_f3) {
        ChannelDef gch;
        gch.pid_combinations = {{21}};
        gch.factors = {1.0};
        channels.push_back(std::move(gch));
    }

    return channels;
}

GridDef load_grid_def(const std::string& path) {
    YAML::Node config = YAML::LoadFile(path);
    GridDef def;

    // Process type
    std::string proc = config["Process"].as<std::string>();
    if (proc == "DIS")        def.process = ProcessType::DIS;
    else if (proc == "SIDIS") def.process = ProcessType::SIDIS;
    else if (proc == "SIA")   def.process = ProcessType::SIA;
    else throw std::runtime_error("Unknown process type: " + proc);

    // Observable (optional, defaults to F2)
    if (config["Observable"]) {
        std::string obs = config["Observable"].as<std::string>();
        if (obs == "F2")       def.observable = Observable::F2;
        else if (obs == "FL")  def.observable = Observable::FL;
        else if (obs == "F3")  def.observable = Observable::F3;
        else throw std::runtime_error("Unknown observable: " + obs);
    }

    // Current (optional, defaults to NC)
    if (config["Current"]) {
        std::string cur = config["Current"].as<std::string>();
        if (cur == "NC")       def.current = Current::NC;
        else if (cur == "CC")  def.current = Current::CC;
        else throw std::runtime_error("Unknown current: " + cur);
    }

    // PID basis
    std::string basis = config["PidBasis"].as<std::string>();
    if (basis == "PDG")       def.pid_basis = PINEAPPL_PID_BASIS_PDG;
    else if (basis == "EVOL") def.pid_basis = PINEAPPL_PID_BASIS_EVOL;
    else throw std::runtime_error("Unknown PID basis: " + basis);

    // Hadron PIDs
    def.hadron_pids = config["HadronPids"].as<std::vector<int>>();

    // Convolution types
    for (const auto& ct : config["ConvolutionTypes"]) {
        std::string s = ct.as<std::string>();
        if (s == "UNPOL_PDF")      def.convolution_types.push_back(PINEAPPL_CONV_TYPE_UNPOL_PDF);
        else if (s == "POL_PDF")   def.convolution_types.push_back(PINEAPPL_CONV_TYPE_POL_PDF);
        else if (s == "UNPOL_FF")  def.convolution_types.push_back(PINEAPPL_CONV_TYPE_UNPOL_FF);
        else if (s == "POL_FF")    def.convolution_types.push_back(PINEAPPL_CONV_TYPE_POL_FF);
        else throw std::runtime_error("Unknown convolution type: " + s);
    }

    // Orders (5-element arrays)
    for (const auto& order : config["Orders"]) {
        auto v = order.as<std::vector<int>>();
        if (v.size() != 5)
            throw std::runtime_error("Each order must have exactly 5 elements");
        def.orders.push_back({static_cast<uint8_t>(v[0]),
                              static_cast<uint8_t>(v[1]),
                              static_cast<uint8_t>(v[2]),
                              static_cast<uint8_t>(v[3]),
                              static_cast<uint8_t>(v[4])});
    }

    // Channels (optional — will be auto-derived in build_grid if absent)
    if (config["Channels"]) {
        for (const auto& ch : config["Channels"]) {
            ChannelDef cdef;
            for (const auto& pids : ch["pids"]) {
                cdef.pid_combinations.push_back(pids.as<std::vector<int>>());
            }
            cdef.factors = ch["factors"].as<std::vector<double>>();
            def.channels.push_back(std::move(cdef));
        }
    }

    // Bins
    for (const auto& bin : config["Bins"]) {
        BinDef bdef;
        bdef.lower = bin["lower"].as<std::vector<double>>();
        bdef.upper = bin["upper"].as<std::vector<double>>();
        def.bins.push_back(std::move(bdef));
    }

    // Normalizations
    def.normalizations = config["Normalizations"].as<std::vector<double>>();

    return def;
}

pineappl_grid* create_grid(const GridDef& def) {
    const std::size_t nb_convolutions = def.convolution_types.size();
    const std::size_t n_bins = def.bins.size();

    // 1. Build channels
    pineappl_channels* channels = pineappl_channels_new(nb_convolutions);
    for (const auto& ch : def.channels) {
        std::size_t combinations = ch.pid_combinations.size();
        std::vector<int32_t> flat_pids;
        for (const auto& combo : ch.pid_combinations) {
            for (int pid : combo) {
                flat_pids.push_back(static_cast<int32_t>(pid));
            }
        }
        pineappl_channels_add(channels, combinations,
                              flat_pids.data(), ch.factors.data());
    }

    // 2. Flatten order params (5 x uint8_t per order)
    std::vector<uint8_t> order_params;
    for (const auto& o : def.orders) {
        order_params.push_back(o.alpha_s);
        order_params.push_back(o.alpha);
        order_params.push_back(o.log_xir);
        order_params.push_back(o.log_xif);
        order_params.push_back(o.log_xia);
    }

    // 3. Placeholder 1D bins [0, 1, 2, ..., n_bins]
    std::vector<double> bin_limits(n_bins + 1);
    std::iota(bin_limits.begin(), bin_limits.end(), 0.0);

    // 4. Build convolution array
    std::vector<pineappl_conv> convolutions(nb_convolutions);
    for (std::size_t i = 0; i < nb_convolutions; i++) {
        convolutions[i].conv_type = def.convolution_types[i];
        convolutions[i].pid = static_cast<int32_t>(def.hadron_pids[i]);
    }

    // 5. Build kinematics based on process type
    std::vector<pineappl_kinematics> kinematics;
    {
        pineappl_kinematics k = {};
        k.tag = PINEAPPL_KINEMATICS_SCALE;
        k.scale = 0;
        kinematics.push_back(k);
    }
    {
        pineappl_kinematics k = {};
        k.tag = PINEAPPL_KINEMATICS_X;
        k.x = 0;
        kinematics.push_back(k);
    }
    if (def.process == ProcessType::SIDIS) {
        pineappl_kinematics k = {};
        k.tag = PINEAPPL_KINEMATICS_X;
        k.x = 1;
        kinematics.push_back(k);
    }

    // 6. Build interpolation specs
    std::vector<pineappl_interp> interps;
    {
        // Scale interpolation
        pineappl_interp si = {};
        si.min = 1e2;
        si.max = 1e8;
        si.nodes = 40;
        si.order = 3;
        si.reweight = PINEAPPL_REWEIGHT_METH_NO_REWEIGHT;
        si.map = PINEAPPL_MAP_APPL_GRID_H0;
        si.interp_meth = PINEAPPL_INTERP_METH_LAGRANGE;
        interps.push_back(si);
    }
    {
        // Momentum fraction interpolation (X)
        pineappl_interp xi = {};
        xi.min = 2e-7;
        xi.max = 1.0;
        xi.nodes = 50;
        xi.order = 3;
        xi.reweight = PINEAPPL_REWEIGHT_METH_APPL_GRID_X;
        xi.map = PINEAPPL_MAP_APPL_GRID_F2;
        xi.interp_meth = PINEAPPL_INTERP_METH_LAGRANGE;
        interps.push_back(xi);

        if (def.process == ProcessType::SIDIS) {
            interps.push_back(xi);  // same spec for X(1)
        }
    }

    // 7. Build scales: {ren=Scale(0), fac=Scale(0), frg=Scale(0)|NoScale}
    pineappl_scale_func_form scales[3] = {};
    // ren = Scale(0)
    scales[0].tag = PINEAPPL_SCALE_FUNC_FORM_SCALE;
    scales[0].scale = 0;
    // fac = Scale(0)
    scales[1].tag = PINEAPPL_SCALE_FUNC_FORM_SCALE;
    scales[1].scale = 0;
    // frg: Scale(0) for SIDIS and SIA, NoScale for DIS
    if (def.process == ProcessType::DIS) {
        scales[2].tag = PINEAPPL_SCALE_FUNC_FORM_NO_SCALE;
    } else {
        scales[2].tag = PINEAPPL_SCALE_FUNC_FORM_SCALE;
        scales[2].scale = 0;
    }

    // 8. Create the grid
    pineappl_grid* grid = pineappl_grid_new2(
        n_bins,
        bin_limits.data(),
        def.orders.size(),
        order_params.data(),
        channels,
        def.pid_basis,
        convolutions.data(),
        interps.size(),
        interps.data(),
        kinematics.data(),
        scales
    );

    // 9. Remap to multi-dimensional bins via pineappl_grid_set_bwfl
    if (!def.bins.empty()) {
        std::size_t bin_dim = def.bins[0].lower.size();
        std::vector<double> bins_lower;
        std::vector<double> bins_upper;
        std::vector<double> norms = def.normalizations;

        // bin-major layout: [dim0_bin0, dim1_bin0, dim0_bin1, dim1_bin1, ...]
        for (const auto& bin : def.bins) {
            for (std::size_t d = 0; d < bin_dim; d++) {
                bins_lower.push_back(bin.lower[d]);
                bins_upper.push_back(bin.upper[d]);
            }
        }

        pineappl_grid_set_bwfl(grid, bins_lower.data(), bins_upper.data(),
                               n_bins, bin_dim, norms.data());
    }

    // 10. Cleanup channels
    pineappl_channels_delete(channels);

    return grid;
}

void set_subgrid(pineappl_grid* grid, std::size_t bin, std::size_t order,
                 std::size_t channel, std::vector<double>& node_values,
                 std::vector<double>& subgrid_array, std::vector<std::size_t>& shape) {
    pineappl_grid_set_subgrid(grid, bin, order, channel,
                              node_values.data(), subgrid_array.data(),
                              shape.data());
}

} // namespace pineapfel
