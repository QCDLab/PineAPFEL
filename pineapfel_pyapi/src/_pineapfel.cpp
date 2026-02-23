#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <pineapfel.h>
#include <pineappl_capi.h>

#include <string>
#include <vector>

namespace py = pybind11;

// RAII wrapper around pineappl_grid*
class PyGrid {
    pineappl_grid *ptr_;

  public:
    explicit PyGrid(pineappl_grid *p) : ptr_(p) {}
    ~PyGrid() {
        if (ptr_) pineappl_grid_delete(ptr_);
    }

    PyGrid(const PyGrid &)            = delete;
    PyGrid &operator=(const PyGrid &) = delete;

    PyGrid(PyGrid &&o) noexcept : ptr_(o.ptr_) { o.ptr_ = nullptr; }

    pineappl_grid *get() const { return ptr_; }

    void           write(const std::string &path) const {
        pineappl_grid_write(ptr_, path.c_str());
    }

    static PyGrid read(const std::string &path) {
        return PyGrid(pineappl_grid_read(path.c_str()));
    }
};

PYBIND11_MODULE(_pineapfel, m) {
    m.doc() = "Python bindings for PineAPFEL";

    py::enum_<pineapfel::ProcessType>(m, "ProcessType")
        .value("DIS", pineapfel::ProcessType::DIS)
        .value("SIDIS", pineapfel::ProcessType::SIDIS)
        .value("SIA", pineapfel::ProcessType::SIA);

    py::enum_<pineapfel::Observable>(m, "Observable")
        .value("F2", pineapfel::Observable::F2)
        .value("FL", pineapfel::Observable::FL)
        .value("F3", pineapfel::Observable::F3);

    py::enum_<pineapfel::Current>(m, "Current")
        .value("NC", pineapfel::Current::NC)
        .value("CC", pineapfel::Current::CC);

    py::enum_<pineapfel::CCSign>(m, "CCSign")
        .value("Plus", pineapfel::CCSign::Plus)
        .value("Minus", pineapfel::CCSign::Minus);

    py::enum_<pineapfel::MassScheme>(m, "MassScheme")
        .value("ZM", pineapfel::MassScheme::ZM)
        .value("FFN", pineapfel::MassScheme::FFN)
        .value("MassiveZero", pineapfel::MassScheme::MassiveZero)
        .value("FONLL", pineapfel::MassScheme::FONLL);

    py::enum_<pineappl_pid_basis>(m, "PidBasis")
        .value("PDG", PINEAPPL_PID_BASIS_PDG)
        .value("EVOL", PINEAPPL_PID_BASIS_EVOL);

    py::enum_<pineappl_conv_type>(m, "ConvolutionType")
        .value("UNPOL_PDF", PINEAPPL_CONV_TYPE_UNPOL_PDF)
        .value("POL_PDF", PINEAPPL_CONV_TYPE_POL_PDF)
        .value("UNPOL_FF", PINEAPPL_CONV_TYPE_UNPOL_FF)
        .value("POL_FF", PINEAPPL_CONV_TYPE_POL_FF);

    py::class_<pineapfel::OrderDef>(m, "OrderDef")
        .def(py::init<>())
        .def(py::init([](int a_s, int a, int r, int f, int ax) {
            return pineapfel::OrderDef{static_cast<uint8_t>(a_s),
                static_cast<uint8_t>(a),
                static_cast<uint8_t>(r),
                static_cast<uint8_t>(f),
                static_cast<uint8_t>(ax)};
        }),
            py::arg("a_s"),
            py::arg("a"),
            py::arg("r"),
            py::arg("f"),
            py::arg("a"))
        .def_readwrite("alpha_s", &pineapfel::OrderDef::alpha_s)
        .def_readwrite("alpha", &pineapfel::OrderDef::alpha)
        .def_readwrite("log_xir", &pineapfel::OrderDef::log_xir)
        .def_readwrite("log_xif", &pineapfel::OrderDef::log_xif)
        .def_readwrite("log_xia", &pineapfel::OrderDef::log_xia);

    py::class_<pineapfel::BinDef>(m, "BinDef")
        .def(py::init<>())
        .def(py::init<std::vector<double>, std::vector<double>>(),
            py::arg("lower"),
            py::arg("upper"))
        .def_readwrite("lower", &pineapfel::BinDef::lower)
        .def_readwrite("upper", &pineapfel::BinDef::upper);

    py::class_<pineapfel::ChannelDef>(m, "ChannelDef")
        .def(py::init<>())
        .def(py::init<std::vector<std::vector<int>>, std::vector<double>>(),
            py::arg("pid_combinations"),
            py::arg("factors"))
        .def_readwrite("pid_combinations",
            &pineapfel::ChannelDef::pid_combinations)
        .def_readwrite("factors", &pineapfel::ChannelDef::factors);

    py::class_<pineapfel::SubGridDef>(m, "SubGridDef")
        .def(py::init<>())
        .def(py::init<int, double, int>(),
            py::arg("n_knots"),
            py::arg("x_min"),
            py::arg("poly_degree"))
        .def_readwrite("n_knots", &pineapfel::SubGridDef::n_knots)
        .def_readwrite("x_min", &pineapfel::SubGridDef::x_min)
        .def_readwrite("poly_degree", &pineapfel::SubGridDef::poly_degree);

    py::class_<pineapfel::TabulationParams>(m, "TabulationParams")
        .def(py::init<>())
        .def(py::init<int, double, int, int>(),
            py::arg("n_points"),
            py::arg("q_min"),
            py::arg("n_steps"),
            py::arg("interp_degree"))
        .def_readwrite("n_points", &pineapfel::TabulationParams::n_points)
        .def_readwrite("q_min", &pineapfel::TabulationParams::q_min)
        .def_readwrite("n_steps", &pineapfel::TabulationParams::n_steps)
        .def_readwrite("interp_degree",
            &pineapfel::TabulationParams::interp_degree);

    py::class_<pineapfel::TheoryCard>(m, "TheoryCard")
        .def(py::init<>())
        .def_readwrite("mu0", &pineapfel::TheoryCard::mu0)
        .def_readwrite("pert_ord", &pineapfel::TheoryCard::pert_ord)
        .def_readwrite("q_ref", &pineapfel::TheoryCard::q_ref)
        .def_readwrite("alpha_qcd_ref", &pineapfel::TheoryCard::alpha_qcd_ref)
        .def_readwrite("quark_thresholds",
            &pineapfel::TheoryCard::quark_thresholds)
        .def_readwrite("flavors", &pineapfel::TheoryCard::flavors)
        .def_readwrite("ckm", &pineapfel::TheoryCard::ckm)
        .def_readwrite("qed", &pineapfel::TheoryCard::qed)
        .def_readwrite("alpha_qed_ref", &pineapfel::TheoryCard::alpha_qed_ref)
        .def_readwrite("lepton_thresholds",
            &pineapfel::TheoryCard::lepton_thresholds)
        .def_readwrite("heavy_quark_masses",
            &pineapfel::TheoryCard::heavy_quark_masses)
        .def_readwrite("mass_nxi", &pineapfel::TheoryCard::mass_nxi)
        .def_readwrite("mass_ximin", &pineapfel::TheoryCard::mass_ximin)
        .def_readwrite("mass_ximax", &pineapfel::TheoryCard::mass_ximax)
        .def_readwrite("mass_intdeg", &pineapfel::TheoryCard::mass_intdeg)
        .def_readwrite("mass_lambda", &pineapfel::TheoryCard::mass_lambda)
        .def_readwrite("mass_imod", &pineapfel::TheoryCard::mass_imod);

    py::class_<pineapfel::OperatorCard>(m, "OperatorCard")
        .def(py::init<>())
        .def_readwrite("xgrid", &pineapfel::OperatorCard::xgrid)
        .def_readwrite("tabulation", &pineapfel::OperatorCard::tabulation)
        .def_readwrite("xi", &pineapfel::OperatorCard::xi);

    py::class_<pineapfel::GridDef>(m, "GridDef")
        .def(py::init<>())
        .def_readwrite("process", &pineapfel::GridDef::process)
        .def_readwrite("observable", &pineapfel::GridDef::observable)
        .def_readwrite("current", &pineapfel::GridDef::current)
        .def_readwrite("cc_sign", &pineapfel::GridDef::cc_sign)
        .def_readwrite("mass_scheme", &pineapfel::GridDef::mass_scheme)
        .def_readwrite("pid_basis", &pineapfel::GridDef::pid_basis)
        .def_readwrite("hadron_pids", &pineapfel::GridDef::hadron_pids)
        .def_readwrite("convolution_types",
            &pineapfel::GridDef::convolution_types)
        .def_readwrite("orders", &pineapfel::GridDef::orders)
        .def_readwrite("channels", &pineapfel::GridDef::channels)
        .def_readwrite("bins", &pineapfel::GridDef::bins)
        .def_readwrite("normalizations", &pineapfel::GridDef::normalizations);

    py::class_<PyGrid>(m, "Grid")
        .def_static("read", &PyGrid::read, py::arg("path"))
        .def("write", &PyGrid::write, py::arg("path"));

    m.def("load_theory_card", &pineapfel::load_theory_card, py::arg("path"));
    m.def("load_operator_card",
        &pineapfel::load_operator_card,
        py::arg("path"));
    m.def("load_grid_def", &pineapfel::load_grid_def, py::arg("path"));

    m.def("derive_channels",
        &pineapfel::derive_channels,
        py::arg("process"),
        py::arg("observable"),
        py::arg("current"),
        py::arg("cc_sign"),
        py::arg("nf_max"));

    m.def(
        "build_grid",
        [](const pineapfel::GridDef       &gd,
            const pineapfel::TheoryCard   &tc,
            const pineapfel::OperatorCard &oc) {
            py::gil_scoped_release rel;
            return PyGrid(pineapfel::build_grid(gd, tc, oc));
        },
        py::arg("grid_def"),
        py::arg("theory"),
        py::arg("op_card"));

    m.def(
        "evolve",
        [](const PyGrid                   &g,
            const pineapfel::TheoryCard   &tc,
            const pineapfel::OperatorCard &oc) {
            py::gil_scoped_release rel;
            return PyGrid(pineapfel::evolve(g.get(), tc, oc));
        },
        py::arg("grid"),
        py::arg("theory"),
        py::arg("op_card"));
}
