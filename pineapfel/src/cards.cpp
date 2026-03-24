#include <pineapfel/cards.h>
#include <stdexcept>
#include <yaml-cpp/yaml.h>

namespace pineapfel {

TheoryCard load_theory_card(const std::string &path) {
    YAML::Node config = YAML::LoadFile(path);
    TheoryCard tc;

    tc.mu0              = config["mu0"].as<double>();
    tc.pert_ord         = config["PerturbativeOrder"].as<int>();
    tc.q_ref            = config["QRef"].as<double>();
    tc.alpha_qcd_ref    = config["AlphaQCDRef"].as<double>();
    tc.quark_thresholds = config["QuarkThresholds"].as<std::vector<double>>();
    tc.flavors          = config["Flavors"].as<std::vector<int>>();
    tc.ckm = config["CKM"] ? config["CKM"].as<std::vector<double>>()
                           : std::vector<double>{0.97428 * 0.97428,
                                 0.22530 * 0.22530,
                                 0.003470 * 0.003470,
                                 0.22520 * 0.22520,
                                 0.97345 * 0.97345,
                                 0.04100 * 0.04100,
                                 0.00862 * 0.00862,
                                 0.04030 * 0.04030,
                                 0.999152 * 0.999152};
    tc.qed = config["QED"] ? config["QED"].as<bool>() : false;
    tc.alpha_qed_ref =
        config["AlphaQEDRef"] ? config["AlphaQEDRef"].as<double>() : 0.0;
    tc.lepton_thresholds =
        config["LeptonThresholds"]
            ? config["LeptonThresholds"].as<std::vector<double>>()
            : std::vector<double>{};

    // Heavy quark masses (default to quark_thresholds padded to 6 if absent)
    // APFEL++ massive initializers require exactly 6 ordered masses.
    if (config["HeavyQuarkMasses"])
        tc.heavy_quark_masses =
            config["HeavyQuarkMasses"].as<std::vector<double>>();
    else {
        tc.heavy_quark_masses = tc.quark_thresholds;
        while (tc.heavy_quark_masses.size() < 6)
            tc.heavy_quark_masses.push_back(172.0);
    }

    // Massive scheme tabulation parameters (all optional)
    if (config["MassNxi"]) tc.mass_nxi = config["MassNxi"].as<int>();
    if (config["MassXiMin"]) tc.mass_ximin = config["MassXiMin"].as<double>();
    if (config["MassXiMax"]) tc.mass_ximax = config["MassXiMax"].as<double>();
    if (config["MassIntDeg"]) tc.mass_intdeg = config["MassIntDeg"].as<int>();
    if (config["MassLambda"])
        tc.mass_lambda = config["MassLambda"].as<double>();
    if (config["MassIMod"]) tc.mass_imod = config["MassIMod"].as<int>();

    return tc;
}

OperatorCard load_operator_card(const std::string &path) {
    YAML::Node   config = YAML::LoadFile(path);
    OperatorCard oc;

    for (const auto &sg : config["xgrid"]) {
        SubGridDef def;
        def.n_knots     = sg["n_knots"].as<int>();
        def.x_min       = sg["x_min"].as<double>();
        def.poly_degree = sg["poly_degree"].as<int>();
        oc.xgrid.push_back(def);
    }
    if (config["zgrid"]) {
        for (const auto &sg : config["zgrid"]) {
            SubGridDef def;
            def.n_knots     = sg["n_knots"].as<int>();
            def.x_min       = sg["x_min"].as<double>();
            def.poly_degree = sg["poly_degree"].as<int>();
            oc.zgrid.push_back(def);
        }
    }

    auto tab                    = config["tabulation"];
    oc.tabulation.n_points      = tab["n_points"].as<int>();
    oc.tabulation.q_min         = tab["q_min"].as<double>();
    oc.tabulation.n_steps       = tab["n_steps"].as<int>();
    oc.tabulation.interp_degree = tab["interp_degree"].as<int>();

    oc.xi                       = config["xi"].as<std::vector<double>>();

    if (config["sidis_mode"]) {
        oc.sidis_mode = config["sidis_mode"].as<std::string>();
        if (oc.sidis_mode != "legacy" && oc.sidis_mode != "bsf_exact") {
            throw std::runtime_error("Unknown sidis_mode: " + oc.sidis_mode);
        }
    }

    if (config["sidis_int_eps"])
        oc.sidis_int_eps = config["sidis_int_eps"].as<double>();
    if (config["sidis_q2_n_intermediate"])
        oc.sidis_q2_n_intermediate =
            config["sidis_q2_n_intermediate"].as<int>();
    if (config["sidis_q2_include_thresholds"])
        oc.sidis_q2_include_thresholds =
            config["sidis_q2_include_thresholds"].as<bool>();
    if (config["sidis_q2_use_bin_centers_only"])
        oc.sidis_q2_use_bin_centers_only =
            config["sidis_q2_use_bin_centers_only"].as<bool>();
    if (config["sidis_z_quad_subdivisions"])
        oc.sidis_z_quad_subdivisions =
            config["sidis_z_quad_subdivisions"].as<int>();

    if (oc.sidis_q2_n_intermediate < 0)
        throw std::runtime_error("sidis_q2_n_intermediate must be >= 0");
    if (oc.sidis_z_quad_subdivisions < 1)
        throw std::runtime_error("sidis_z_quad_subdivisions must be >= 1");

    return oc;
}

} // namespace pineapfel
