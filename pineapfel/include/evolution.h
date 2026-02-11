#pragma once

#include <cards.h>
#include <pineappl_capi.h>

namespace pineapfel {

// Evolve a PineAPPL grid into an FK table.
pineappl_grid* evolve(pineappl_grid* grid, const TheoryCard& theory, const OperatorCard& op_card);

} // namespace pineapfel
