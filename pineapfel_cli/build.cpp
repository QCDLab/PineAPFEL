#include <NeoPDF.hpp>
#include <apfel/apfelxx.h>
#include <pineapfel/pineapfel.h>
#include <pineappl_capi.h>

#include <iostream>
#include <string>

void print_usage(const char *prog) {
    std::cerr
        << "Usage: " << prog
        << " <grid.yaml> <theory.yaml> <operator.yaml> [-o output.pineappl.lz4]"
        << std::endl;
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        print_usage(argv[0]);
        return 1;
    }

    std::string grid_path   = argv[1];
    std::string theory_path = argv[2];
    std::string op_path     = argv[3];
    std::string output_path;

    for (int i = 4; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-o" && i + 1 < argc) {
            output_path = argv[++i];
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            print_usage(argv[0]);
            return 1;
        }
    }

    if (output_path.empty()) {
        // Default: replace .yaml with .pineappl.lz4
        auto pos = grid_path.rfind(".yaml");
        if (pos != std::string::npos) {
            output_path = grid_path.substr(0, pos) + ".pineappl.lz4";
        } else {
            output_path = grid_path + ".pineappl.lz4";
        }
    }

    std::cout << "Grid card: " << grid_path << std::endl;
    std::cout << "Theory:    " << theory_path << std::endl;
    std::cout << "Operator:  " << op_path << std::endl;
    std::cout << "Output:    " << output_path << std::endl;

    apfel::Timer   t;

    auto           grid_def = pineapfel::load_grid_def(grid_path);
    auto           theory   = pineapfel::load_theory_card(theory_path);
    auto           op_card  = pineapfel::load_operator_card(op_path);

    pineappl_grid *grid     = pineapfel::build_grid(grid_def, theory, op_card);

    // Apply FixConvolutions entries in order, each folding in one PDF/FF set.
    pineappl_grid *result   = grid;
    for (const auto &fc : grid_def.fix_convolutions) {
        std::cout << "Fixing convolution " << fc.index << " with " << fc.pdf_set
                  << " (member " << fc.pdf_member << ", xi=" << fc.xi << ")"
                  << std::endl;

        neopdf::NeoPDF pdf{fc.pdf_set, static_cast<std::size_t>(fc.pdf_member)};

        auto cb = [](int32_t pid, double x, double q2, void *s) -> double {
            return static_cast<neopdf::NeoPDF *>(s)->xfxQ2(pid, x, q2);
        };
        pineappl_grid *fixed = pineappl_grid_fix_convolution(result,
            fc.index,
            cb,
            static_cast<void *>(&pdf),
            fc.xi);

        pineappl_grid_delete(result);
        result = fixed;
    }

    pineappl_grid_optimize(result);
    pineappl_grid_write(result, output_path.c_str());
    std::cout << "Grid written to: " << output_path << std::endl;

    pineappl_grid_delete(result);
    t.stop();

    return 0;
}
