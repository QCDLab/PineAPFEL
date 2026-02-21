#pragma once

#include <cards.h>
#include <grid_gen.h>
#include <pineappl_capi.h>

namespace pineapfel {

// Build and fill a PineAPPL grid with APFEL++ coefficient functions.
// Supports DIS, SIA, and SIDIS processes with NC and CC currents in
// the zero-mass scheme.
pineappl_grid *build_grid(const GridDef &grid_def,
    const TheoryCard                    &theory,
    const OperatorCard                  &op_card);

} // namespace pineapfel
