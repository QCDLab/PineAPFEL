#pragma once

#include <apfel/apfelxx.h>
#include <apfel/sidisbuilder.h>
#include <functional>
#include <vector>

// Compute SIDIS F2 reference values using APFEL++ InitializeSidisObjects (exact
// NNLO). Returns per-bin reference values for comparison with PineAPPL
// convolution.
// sobj must be pre-initialised (from init_sidis_nnlo) and built on the same
// grid g used when filling the PineAPPL grid.
std::vector<double> compute_sidis_reference(const apfel::Grid &g,
    const apfel::SidisNNLOObjects                             &sobj,
    const std::vector<double>                                 &thresholds,
    const std::vector<double>                                 &q2_nodes,
    const std::vector<std::vector<double>> &bin_x_bounds,  // {lo, hi} per bin
    const std::vector<std::vector<double>> &bin_z_bounds,  // {lo, hi} per bin
    const std::vector<double>              &charges_dummy, // unused placeholder
    int                                     max_alpha_s,
    std::function<double(int, double)>      toy_f,
    std::function<double(double)>           alphas_func);

// Compute polarized SIDIS G1 reference values using APFEL++
// InitializeSidisObjects. Returns per-bin reference values for comparison with
// PineAPPL convolution.
// sobj must be pre-initialised (from init_sidis_nnlo) and built on the same
// grid g used when filling the PineAPPL grid.
std::vector<double> compute_sidis_pol_reference(const apfel::Grid &g,
    const apfel::SidisNNLOObjects                                 &sobj,
    const std::vector<double>                                     &thresholds,
    const std::vector<double>                                     &q2_nodes,
    const std::vector<std::vector<double>> &bin_x_bounds, // {lo, hi} per bin
    const std::vector<std::vector<double>> &bin_z_bounds, // {lo, hi} per bin
    int                                     max_alpha_s,
    std::function<double(int, double)>      toy_f,
    std::function<double(double)>           alphas_func);

// Compute SIDIS F2 NNLO-only (alpha_s^2) reference values by directly
// evaluating each DoubleOperator via dop.MultiplyFirstBy(pdf).Evaluate(x)*ff.
// This is the gold-standard comparison for the eval_double_op_column kernel.
// sobj must be pre-initialised (from init_sidis_nnlo); nf_max must match
// derive_channels.
std::vector<double> compute_sidis_nnlo_reference(const apfel::Grid &g,
    const apfel::SidisNNLOObjects                                  &sobj,
    const std::vector<double>                                      &thresholds,
    int                                                             nf_max,
    const std::vector<double>                                      &q2_nodes,
    const std::vector<std::vector<double>> &bin_x_bounds,
    const std::vector<std::vector<double>> &bin_z_bounds,
    std::function<double(int, double)>      toy_f,
    std::function<double(double)>           alphas_func);
