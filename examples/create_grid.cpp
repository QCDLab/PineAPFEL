#include <pineapfel.h>
#include <pineappl_capi.h>

#include <iostream>

int main() {
    // Load cards
    auto grid_def = pineapfel::load_grid_def("runcards/grid_dis.yaml");
    auto theory   = pineapfel::load_theory_card("runcards/theory.yaml");
    auto op_card  = pineapfel::load_operator_card("runcards/operator.yaml");

    std::cout << "Grid definition:" << std::endl;
    std::cout << "  Process:    "
              << (grid_def.process == pineapfel::ProcessType::DIS ? "DIS"
                                                                  : "SIA")
              << std::endl;
    std::cout << "  Bins:       " << grid_def.bins.size() << std::endl;
    std::cout << "  Orders:     " << grid_def.orders.size() << std::endl;
    std::cout << "  Channels:   " << grid_def.channels.size() << std::endl;

    // Build and fill the grid with APFEL++ coefficient functions
    auto       *grid   = pineapfel::build_grid(grid_def, theory, op_card);

    // Write the filled grid to disk
    const char *output = "dis_test.pineappl.lz4";
    pineappl_grid_write(grid, output);
    std::cout << "Grid written to: " << output << std::endl;

    pineappl_grid_delete(grid);
    return 0;
}
