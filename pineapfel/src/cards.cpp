#include <cards.h>
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
    tc.ckm              = config["CKM"]
                              ? config["CKM"].as<std::vector<double>>()
                              : std::vector<double>{
                                    0.97428 * 0.97428, 0.22530 * 0.22530,
                                    0.003470 * 0.003470, 0.22520 * 0.22520,
                                    0.97345 * 0.97345, 0.04100 * 0.04100,
                                    0.00862 * 0.00862, 0.04030 * 0.04030,
                                    0.999152 * 0.999152};
    tc.qed              = config["QED"] ? config["QED"].as<bool>() : false;
    tc.alpha_qed_ref =
        config["AlphaQEDRef"] ? config["AlphaQEDRef"].as<double>() : 0.0;
    tc.lepton_thresholds =
        config["LeptonThresholds"]
            ? config["LeptonThresholds"].as<std::vector<double>>()
            : std::vector<double>{};

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

    auto tab                    = config["tabulation"];
    oc.tabulation.n_points      = tab["n_points"].as<int>();
    oc.tabulation.q_min         = tab["q_min"].as<double>();
    oc.tabulation.n_steps       = tab["n_steps"].as<int>();
    oc.tabulation.interp_degree = tab["interp_degree"].as<int>();

    oc.xi                       = config["xi"].as<std::vector<double>>();

    return oc;
}

} // namespace pineapfel
