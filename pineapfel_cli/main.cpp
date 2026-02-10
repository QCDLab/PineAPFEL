#include <pineapfel/pineapfel.h>
#include <pineappl_capi.h>
#include <apfel/apfelxx.h>

#include <iostream>
#include <string>
#include <vector>

void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " <grid.pineappl.lz4> <theory.yaml> <operator.yaml> [-o output.pineappl.lz4]"
              << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        print_usage(argv[0]);
        return 1;
    }

    std::string grid_path   = argv[1];
    std::string theory_path = argv[2];
    std::string op_path     = argv[3];
    std::string output_path;

    // Parse optional -o flag
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

    // Default output path
    if (output_path.empty()) {
        auto pos = grid_path.find(".pineappl.lz4");
        if (pos != std::string::npos) {
            output_path = grid_path.substr(0, pos) + ".fk.pineappl.lz4";
        } else {
            output_path = grid_path + ".fk.pineappl.lz4";
        }
    }

    std::cout << "Grid:     " << grid_path << std::endl;
    std::cout << "Theory:   " << theory_path << std::endl;
    std::cout << "Operator: " << op_path << std::endl;
    std::cout << "Output:   " << output_path << std::endl;

    apfel::Timer t;
    auto theory  = pineapfel::load_theory_card(theory_path);
    auto op_card = pineapfel::load_operator_card(op_path);
    auto* grid = pineappl_grid_read(grid_path.c_str());

    auto* fktable = pineapfel::evolve(grid, theory, op_card);
    pineappl_grid_write(fktable, output_path.c_str());
    std::cout << "FK table written to: " << output_path << std::endl;

    pineappl_grid_delete(grid);
    pineappl_grid_delete(fktable);
    t.stop();

    return 0;
}
