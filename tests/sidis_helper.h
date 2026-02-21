#pragma once

#include <apfel/apfelxx.h>
#include <functional>
#include <vector>

// Compute SIDIS F2 reference values using APFEL++ InitializeSIDIS.
// Returns per-bin reference values for comparison with PineAPPL convolution.
// toy_f(pid, x) returns f(x) (NOT x*f).
std::vector<double> compute_sidis_reference(const apfel::Grid &g,
    const std::vector<double>                                 &thresholds,
    const std::vector<double>                                 &q2_nodes,
    const std::vector<std::vector<double>> &bin_x_bounds,  // {lo, hi} per bin
    const std::vector<std::vector<double>> &bin_z_bounds,  // {lo, hi} per bin
    const std::vector<double>              &charges_dummy, // unused placeholder
    int                                     max_alpha_s,
    std::function<double(int, double)>      toy_f,
    std::function<double(double)>           alphas_func);
