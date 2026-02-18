  Current State

  PineAPFEL already has two well-separated concerns:
  1. Grid creation — create_grid() + set_subgrid() (from YAML or programmatic GridDef)
  2. Evolution — evolve() turns a PineAPPL grid into an FK table using APFEL++ DGLAP operators

  The missing piece is step 1.5: using APFEL++ to compute the actual coefficient functions and fill the subgrids — currently done
  with random numbers in the example.

  ---
  Recommended Architecture

  A single new function per process, driven entirely by runcards

  The cleanest approach is a new module (e.g. fill.h / fill.cpp) that provides:

  // The main user-facing function
  pineappl_grid* build_grid(
      const GridDef&      grid_def,
      const TheoryCard&   theory,
      const OperatorCard& op_card
  );

  Internally, it would dispatch on grid_def.process:

  ┌─────────┬─────────────────────────────────────┬───────────────┬──────────────────────┐
  │ Process │         APFEL++ initializer         │ Subgrid shape │      Kinematics      │
  ├─────────┼─────────────────────────────────────┼───────────────┼──────────────────────┤
  │ DIS     │ InitializeF2NCObjectsZM() etc.      │ (nQ, nx)      │ Scale(0), X(0)       │
  ├─────────┼─────────────────────────────────────┼───────────────┼──────────────────────┤
  │ SIA     │ InitializeF2NCObjectsZMT() etc.     │ (nQ, nz)      │ Scale(0), X(0)       │
  ├─────────┼─────────────────────────────────────┼───────────────┼──────────────────────┤
  │ SIDIS   │ SIDIS DoubleExpression coefficients │ (nQ, nx, nz)  │ Scale(0), X(0), X(1) │
  └─────────┴─────────────────────────────────────┴───────────────┴──────────────────────┘

  Internal flow

  User provides 3 YAML files
         ↓
  load_grid_def("grid_dis.yaml")   → GridDef
  load_theory_card("theory.yaml")  → TheoryCard
  load_operator_card("operator.yaml") → OperatorCard
         ↓
  build_grid(grid_def, theory, op_card)
    1. Build apfel::Grid from op_card.xgrid (same grid used later for evolution)
    2. Get Q² nodes from tabulation params (or from bin edges)
    3. Get x/z nodes from the apfel::Grid joint grid
    4. Dispatch to process-specific filler:
       ├─ DIS:   init StructureFunctionObjects (ZM or massive)
       │         for each (order, channel, bin):
       │           evaluate C(x) ⊗ basis at each Q² node → fill 2D subgrid
       ├─ SIA:   same but with timelike initializers + FF convolution type
       └─ SIDIS: init SIDIS coefficient DoubleExpressions
                for each (order, channel, bin):
                  evaluate C(x,z) at each Q² node → fill 3D subgrid
    5. create_grid(grid_def) + set_subgrid() for all entries
    6. Return filled grid


  What the grid YAML cards need to encode additionally

  The current grid cards define bins, channels, and orders but not which structure functions to compute (F2, FL, F3, g1...). You'd
  want to add something like:

  # grid_dis.yaml
  process: DIS
  observable: F2     # or FL, F3, g1, etc.
  scheme: ZM         # or MASSIVE
  current: NC        # or CC

  This tells the filler which Initialize*Objects*() function to call.

  Why this design

  1. User simplicity — The entire workflow becomes:
  auto grid_def = pineapfel::load_grid_def("grid_dis.yaml");
  auto theory   = pineapfel::load_theory_card("theory.yaml");
  auto op_card  = pineapfel::load_operator_card("operator.yaml");
  auto* grid    = pineapfel::build_grid(grid_def, theory, op_card);
  auto* fk      = pineapfel::evolve(grid, theory, op_card);
  1. Or even a single CLI command: pineapfel-build grid_dis.yaml theory.yaml operator.yaml
  2. Separation stays clean — build_grid fills coefficient functions, evolve does DGLAP. They share the same x-grid (from
  OperatorCard) which is important for consistency.
  3. SIDIS is the hard part — DIS and SIA are structurally identical (1 convolution, 2D subgrids, Expression-based coefficients).
  SIDIS needs DoubleExpression and 3D subgrids but the PineAPFEL grid infrastructure already supports this via ProcessType::SIDIS
  with {Scale(0), X(0), X(1)}.
  4. Channels map naturally — APFEL++'s DISNCBasis / DISCCBasis decompose into non-singlet/singlet/gluon components, which
  correspond directly to the channel definitions in the grid YAML.

  Key implementation detail

  The tricky part is the mapping between PineAPPL channels and APFEL++ basis elements. The grid YAML defines channels by PID
  combinations (e.g., {2, -2} for u+ubar). The APFEL++ coefficient functions are computed in the {CNS, CS, CG} operand basis and
  then rotated via DISNCBasis. You'll need a translation layer that:

  - Takes the APFEL++ StructureFunctionObjects (which contain C0, C1, C2, C3 as Set<Operator>)
  - Evaluates each operator on the x-grid nodes
  - Maps the result to the PID-channel structure expected by PineAPPL

  This is essentially the inverse of what BuildStructureFunctions() does — instead of convolving with actual PDFs, you extract the
  bare operator values on the interpolation grid.


 Ready to code?

 Here is Claude's plan:
╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
 Plan: APFEL++ Coefficient Function Filling for PineAPFEL

 Context

 PineAPFEL currently has two modules: grid creation (grid_gen) and DGLAP evolution (evolution). The missing piece is filling
 PineAPPL grids with APFEL++ coefficient functions so users only provide YAML runcards. Currently examples/create_grid.cpp fills
 subgrids with random data.

 Goal: A new build_grid() function that reads 3 YAML cards (grid, theory, operator) and returns a PineAPPL grid filled with
 analytically computed coefficient functions from APFEL++. Scope: DIS + SIA, zero-mass scheme, PDG basis channels.

 ---
 Step 1: Extend GridDef and YAML cards

 Files: pineapfel/include/grid_gen.h, pineapfel/src/grid_gen.cpp

 Add new fields to GridDef:

 // grid_gen.h — new enums and fields
 enum class Observable { F2, FL, F3 };
 enum class Current    { NC, CC };

 struct GridDef {
     // ... existing fields ...
     Observable observable = Observable::F2;
     Current    current    = Current::NC;
 };

 Update load_grid_def() in grid_gen.cpp to parse new YAML fields:

 # grid_dis.yaml — new fields
 Observable: F2
 Current: NC

 Update runcards/grid_dis.yaml and runcards/grid_sia.yaml with these new fields.

 ---
 Step 2: New module fill.h / fill.cpp

 New files: pineapfel/include/fill.h, pineapfel/src/fill.cpp

 Public API

 // fill.h
 #pragma once
 #include <cards.h>
 #include <grid_gen.h>
 #include <pineappl_capi.h>

 namespace pineapfel {

 pineappl_grid* build_grid(
     const GridDef&      grid_def,
     const TheoryCard&   theory,
     const OperatorCard& op_card
 );

 } // namespace pineapfel

 Internal flow of build_grid()

 build_grid(grid_def, theory, op_card)
 │
 ├─ 1. Build APFEL++ x-grid from op_card.xgrid
 │     (reuse same pattern as evolution.cpp:222-227)
 │
 ├─ 2. Select APFEL++ initializer based on (process, observable, current)
 │     DIS + F2 + NC → InitializeF2NCObjectsZM(g, thresholds)
 │     DIS + FL + NC → InitializeFLNCObjectsZM(g, thresholds)
 │     DIS + F3 + NC → InitializeF3NCObjectsZM(g, thresholds)
 │     SIA + F2 + NC → InitializeF2NCObjectsZMT(g, thresholds)
 │     SIA + FL + NC → InitializeFLNCObjectsZMT(g, thresholds)
 │
 ├─ 3. Create the empty PineAPPL grid via create_grid(grid_def)
 │
 ├─ 4. Determine subgrid nodes:
 │     x_nodes = g.GetJointGrid().GetGrid()  (APFEL++ joint grid)
 │     Q²_nodes = derive from bin edges (union of all bin Q² ranges)
 │
 ├─ 5. For each Q² node:
 │     │  Q = sqrt(Q²)
 │     │  nf = apfel::NF(Q, thresholds)
 │     │  charges = apfel::ElectroWeakCharges(Q, false)  [DIS: spacelike]
 │     │         or apfel::ElectroWeakCharges(Q, true)   [SIA: timelike]
 │     │  FObj = initializer(Q, charges)  → StructureFunctionObjects
 │     │
 │     │  Extract raw operators at each perturbative order n:
 │     │    C_NS^{(n)} = FObj.Cn.at(0).GetObjects().at(DISNCBasis::CNS)
 │     │    C_CS^{(n)} = FObj.Cn.at(0).GetObjects().at(DISNCBasis::CS)
 │     │    C_CG^{(n)} = FObj.Cn.at(0).GetObjects().at(DISNCBasis::CG)
 │     │    C_PS^{(n)} = (C_CS - C_NS) [pure singlet, zero at LO/NLO]
 │     │
 │     └─ 6. For each (bin, order, channel):
 │           │  x_B = bin center in x dimension (geometric mean of lower/upper)
 │           │  order_idx = map OrderDef.alpha_s to APFEL++ order (0→C0, 1→C1, 2→C2)
 │           │
 │           │  Channel mapping (PDG PIDs → operator combination):
 │           │  ┌──────────────────────────────────────────────────────┐
 │           │  │ Channel [[q],[-q]] (quark q + antiquark):           │
 │           │  │   C_channel = e_q² * C_NS + (Σ_i e_i²) * C_PS     │
 │           │  │                                                      │
 │           │  │ Channel [[21]] (gluon):                              │
 │           │  │   C_channel = (Σ_q e_q²) * C_CG                    │
 │           │  └──────────────────────────────────────────────────────┘
 │           │
 │           │  Evaluate: dist = C_channel.Evaluate(x_B)
 │           │  subgrid[Q²_idx * nx + x_idx] = dist.GetDistributionJointGrid()[x_idx]
 │           │
 │           └─ Call set_subgrid(grid, bin, order, channel, ...)
 │
 └─ Return filled grid


 Channel → operator mapping detail

 For DIS/SIA F2 NC at perturbative order n, using Set<Operator> from FObj.Cn.at(k=0):

 // Extract raw operators (charge-independent)
 auto& ops = FObj.C{n}.at(0).GetObjects();
 Operator CNS = ops.at(DISNCBasis::CNS);
 Operator CS  = ops.at(DISNCBasis::CS);
 Operator CG  = ops.at(DISNCBasis::CG);
 Operator CPS = CS - CNS;  // pure singlet (zero at LO, NLO)

 // For quark channel with PID q:
 int q_idx = abs(q) - 1;  // 0-indexed flavor
 double e_q_sq = charges[q_idx];
 double sum_e_sq = accumulate(charges for nf active flavors);
 Operator C_quark = e_q_sq * CNS + sum_e_sq * CPS;

 // For gluon channel:
 Operator C_gluon = sum_e_sq * CG;

 Note on F3: For F3, the non-singlet coefficient is the minus-type (C_NS^-), and there is no singlet/gluon contribution. The
 channel mapping simplifies to C_channel = e_q² * C_NS for quarks.

 Q² node selection

 Derive Q² nodes from bin edges to ensure adequate coverage:
 - Collect all unique Q² values from bin lower/upper edges
 - Add intermediate points (geometric spacing) within each bin's Q² range
 - Include points at quark thresholds if they fall within the range

 Subgrid node_values and shape

 Following the pattern in examples/create_grid.cpp:

 // DIS: node_values = [Q²_0, ..., Q²_{nq-1}, x_0, ..., x_{nx-1}]
 //       shape = {nq, nx}
 // SIA:  node_values = [Q²_0, ..., Q²_{nq-1}, z_0, ..., z_{nz-1}]
 //       shape = {nq, nz}

 The x/z nodes come from g.GetJointGrid().GetGrid() (APFEL++ joint grid, same grid used for evolution). This ensures consistency
 between grid filling and later evolution.

 ---
 Step 3: Update YAML runcards

 runcards/grid_dis.yaml — add:
 Observable: F2
 Current: NC

 runcards/grid_sia.yaml — add:
 Observable: F2
 Current: NC

 ---
 Step 4: Update build system

 meson.build — add 'pineapfel/src/fill.cpp' to pineapfel_lib sources, add fill.h to installed headers.

 CMakeLists.txt — add pineapfel/src/fill.cpp to add_library(pineapfel ...).

 ---
 Step 5: Update pineapfel.h convenience header

 Add #include <fill.h> to pineapfel/include/pineapfel.h.

 ---
 Step 6: Add CLI command for grid building

 File: pineapfel_cli/main.cpp

 Add a build sub-command (or a new executable pineapfel-build):

 ./pineapfel-build <grid.yaml> <theory.yaml> <operator.yaml> [-o output.pineappl.lz4]


 This calls build_grid() and writes the result. Can be kept simple — possibly just extend the existing CLI or add a second
 executable.

 ---
 Step 7: Update example

 File: examples/create_grid.cpp

 Replace the random-fill example with an APFEL++ fill example, demonstrating:
 auto grid_def = pineapfel::load_grid_def("runcards/grid_dis.yaml");
 auto theory   = pineapfel::load_theory_card("runcards/theory.yaml");
 auto op_card  = pineapfel::load_operator_card("runcards/operator.yaml");
 auto* grid    = pineapfel::build_grid(grid_def, theory, op_card);
 pineappl_grid_write(grid, "dis_f2.pineappl.lz4");
 pineappl_grid_delete(grid);

 ---
 Files to create/modify

 ┌───────────────────────────────┬─────────────────────────────────────────────────┐
 │             File              │                     Action                      │
 ├───────────────────────────────┼─────────────────────────────────────────────────┤
 │ pineapfel/include/fill.h      │ Create — build_grid() declaration               │
 ├───────────────────────────────┼─────────────────────────────────────────────────┤
 │ pineapfel/src/fill.cpp        │ Create — Core implementation (~200-300 lines)   │
 ├───────────────────────────────┼─────────────────────────────────────────────────┤
 │ pineapfel/include/grid_gen.h  │ Edit — Add Observable, Current enums to GridDef │
 ├───────────────────────────────┼─────────────────────────────────────────────────┤
 │ pineapfel/src/grid_gen.cpp    │ Edit — Parse new YAML fields in load_grid_def() │
 ├───────────────────────────────┼─────────────────────────────────────────────────┤
 │ pineapfel/include/pineapfel.h │ Edit — Add #include <fill.h>                    │
 ├───────────────────────────────┼─────────────────────────────────────────────────┤
 │ runcards/grid_dis.yaml        │ Edit — Add Observable, Current                  │
 ├───────────────────────────────┼─────────────────────────────────────────────────┤
 │ runcards/grid_sia.yaml        │ Edit — Add Observable, Current                  │
 ├───────────────────────────────┼─────────────────────────────────────────────────┤
 │ meson.build                   │ Edit — Add fill.cpp source and fill.h header    │
 ├───────────────────────────────┼─────────────────────────────────────────────────┤
 │ CMakeLists.txt                │ Edit — Add fill.cpp source                      │
 ├───────────────────────────────┼─────────────────────────────────────────────────┤
 │ pineapfel_cli/main.cpp        │ Edit — Add build sub-command or new executable  │
 ├───────────────────────────────┼─────────────────────────────────────────────────┤
 │ examples/create_grid.cpp      │ Edit — Use build_grid() instead of random data  │
 └───────────────────────────────┴─────────────────────────────────────────────────┘

 ---
 Verification

 1. Build: Compile with both meson and cmake, ensure no errors
 2. Run example: Execute the updated create_grid example — should produce a .pineappl.lz4 grid with non-zero, physically meaningful
  subgrid values
 3. Evolution roundtrip: Run pineapfel-evolve on the generated grid to produce an FK table — should complete without errors
 4. Sanity check: Use pineappl CLI to inspect the grid (pineappl read <grid>) and verify:
   - Correct number of bins, orders, channels
   - Non-zero subgrid entries
   - Reasonable numerical values (coefficient functions are O(1) at LO)
 5. Cross-check with APFEL++: Compare FK table convolution results with direct APFEL++ BuildStructureFunctions() computation at a
 few (x, Q²) points — values should agree within interpolation precision


● All tasks are complete. Here's a summary of what was implemented:

  Summary

  New files created

  - pineapfel/include/fill.h — Declares build_grid(GridDef, TheoryCard, OperatorCard) returning a filled pineappl_grid*
  - pineapfel/src/fill.cpp — Core implementation (~230 lines) that:
    a. Builds an APFEL++ x-space grid from the operator card
    b. Selects the correct initializer (DIS/SIA × F2/FL/F3 × NC)
    c. Creates an empty PineAPPL grid via create_grid()
    d. Derives Q² nodes from bin edges with geometric intermediate points
    e. Precomputes coefficient function operators at each Q² node
    f. For each (bin, order, channel), constructs the per-channel coefficient function operator and fills the subgrid
  - pineapfel_cli/build.cpp — New pineapfel-build CLI executable

  Modified files

  ┌───────────────────────────────┬──────────────────────────────────────────────────────────────────────┐
  │             File              │                                Change                                │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────┤
  │ pineapfel/include/grid_gen.h  │ Added Observable (F2/FL/F3) and Current (NC/CC) enums to GridDef     │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────┤
  │ pineapfel/src/grid_gen.cpp    │ Parse Observable and Current from YAML (optional, defaults to F2/NC) │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────┤
  │ pineapfel/include/pineapfel.h │ Added #include <fill.h>                                              │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────┤
  │ runcards/grid_dis.yaml        │ Added Observable: F2, Current: NC                                    │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────┤
  │ runcards/grid_sia.yaml        │ Added Observable: F2, Current: NC                                    │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────┤
  │ meson.build                   │ Added fill.cpp source, fill.h header, pineapfel-build executable     │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────┤
  │ CMakeLists.txt                │ Added fill.cpp source, pineapfel-build executable                    │
  ├───────────────────────────────┼──────────────────────────────────────────────────────────────────────┤
  │ examples/create_grid.cpp      │ Now uses build_grid() instead of random data                         │
  └───────────────────────────────┴──────────────────────────────────────────────────────────────────────┘

  Key physics formulas

  The coefficient function decomposition from evolution basis to physical (PDG) channels:
  - Quark channel [[q],[-q]]: C_q = e_q² × C_NS + (ΣCh/6) × (C_S − C_NS)
  - Gluon channel [[21]]: C_g = ΣCh × C_G
  - F3 (valence-only): C_q = e_q² × C_NS, C_g = 0

  The factor of 6 matches the hardcoded normalization in APFEL++'s DISNCBasis convolution map (where _rules[SIGMA] = {CS, SIGMA,
  SumCh/6}).
