#include <pineapfel.h>
#include <pineappl_capi.h>

#include <cstddef>
#include <iostream>
#include <random>
#include <vector>

int main() {
    // ---------------------------------------------------------------
    // Option A: Load the grid definition from a YAML file
    // ---------------------------------------------------------------
    //   auto def = pineapfel::load_grid_def("runcards/grid_dis.yaml");

    // ---------------------------------------------------------------
    // Option B: Build the grid definition programmatically
    // ---------------------------------------------------------------
    pineapfel::GridDef def;
    def.process = pineapfel::ProcessType::DIS;
    def.pid_basis = PINEAPPL_PID_BASIS_PDG;
    def.hadron_pids = { 2212 };
    def.convolution_types = { PINEAPPL_CONV_TYPE_UNPOL_PDF };

    // O(alpha_s^2) with no log terms
    def.orders = { {2, 0, 0, 0, 0} };

    // Three channels: u+ubar, d+dbar, gluon
    def.channels = {
        { /*pids=*/{{2}, {-2}},  /*factors=*/{1.0, 1.0} },
        { /*pids=*/{{1}, {-1}},  /*factors=*/{1.0, 1.0} },
        { /*pids=*/{{21}},       /*factors=*/{1.0} },
    };

    // Two bins in (Q^2, x)
    def.bins = {
        { /*lower=*/{10.0, 0.001},  /*upper=*/{100.0, 0.01} },
        { /*lower=*/{100.0, 0.01},  /*upper=*/{1000.0, 0.1} },
    };
    def.normalizations = { 1.0, 1.0 };

    // ---------------------------------------------------------------
    // Create the empty PineAPPL grid
    // ---------------------------------------------------------------
    auto* grid = pineapfel::create_grid(def);

    std::cout << "Created DIS grid" << std::endl;
    std::cout << "  Bins:     " << pineappl_grid_bin_count(grid) << std::endl;
    std::cout << "  Orders:   " << def.orders.size()             << std::endl;
    std::cout << "  Channels: " << def.channels.size()           << std::endl;

    // ---------------------------------------------------------------
    // Fill every (bin, order, channel) with a random subgrid
    // ---------------------------------------------------------------
    //
    // For DIS the kinematics are {Scale(0), X(0)}, so:
    //   node_values = [q2_node_0 ... q2_node_{nq-1}, x_node_0 ... x_node_{nx-1}]
    //   shape       = [nq, nx]
    //   subgrid     = flat array of size nq * nx (row-major)

    // These information will eventually come from APFEL++
    std::vector<double> q2_nodes = { 1e2, 5e2, 1e3, 5e3, 1e4 };
    std::vector<double> x_nodes  = { 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1 };
    const std::size_t nx = x_nodes.size();
    const std::size_t nq = q2_nodes.size();

    // Concatenate: first all scale nodes, then all x nodes
    std::vector<double> node_values;
    node_values.insert(node_values.end(), q2_nodes.begin(), q2_nodes.end());
    node_values.insert(node_values.end(), x_nodes.begin(), x_nodes.end());

    std::vector<std::size_t> shape = {nq, nx};

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (std::size_t bin = 0; bin < def.bins.size(); bin++) {
        for (std::size_t order = 0; order < def.orders.size(); order++) {
            for (std::size_t channel = 0; channel < def.channels.size(); channel++) {
                std::vector<double> subgrid(nq * nx);
                for (auto& v : subgrid) v = dist(rng);

                pineapfel::set_subgrid(
                    grid,
                    bin,
                    order,
                    channel,
                    node_values,
                    subgrid,
                    shape
                );
            }
        }
    }

    // ---------------------------------------------------------------
    // Write the filled grid to disk
    // ---------------------------------------------------------------
    const char* output = "dis_test.pineappl.lz4";
    pineappl_grid_write(grid, output);
    std::cout << "Grid written to: " << output << std::endl;

    pineappl_grid_delete(grid);
    return 0;
}
