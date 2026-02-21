#pragma once

#include <pineappl_capi.h>

#include <cstdint>
#include <string>
#include <vector>

namespace pineapfel {

enum class ProcessType { DIS, SIDIS, SIA };
enum class Observable { F2, FL, F3 };
enum class Current { NC, CC };
enum class CCSign { Plus, Minus };

struct OrderDef {
    uint8_t alpha_s, alpha, log_xir, log_xif, log_xia;
};

struct ChannelDef {
    std::vector<std::vector<int>>
                        pid_combinations; // nb_convolutions PIDs per combo
    std::vector<double> factors;          // one per combination
};

struct BinDef {
    std::vector<double> lower; // one value per dimension
    std::vector<double> upper;
};

struct GridDef {
    ProcessType                     process;
    Observable                      observable = Observable::F2;
    Current                         current    = Current::NC;
    CCSign                          cc_sign    = CCSign::Plus;
    pineappl_pid_basis              pid_basis;
    std::vector<int>                hadron_pids;
    std::vector<pineappl_conv_type> convolution_types;
    std::vector<OrderDef>           orders;
    std::vector<ChannelDef>         channels;
    std::vector<BinDef>             bins;           // one BinDef per bin
    std::vector<double>             normalizations; // one per bin
};

GridDef                 load_grid_def(const std::string &path);
std::vector<ChannelDef> derive_channels(ProcessType process,
    Observable                                      observable,
    Current                                         current,
    CCSign                                          cc_sign,
    int                                             nf_max);
pineappl_grid          *create_grid(const GridDef &def);
void                    set_subgrid(pineappl_grid *grid,
                       std::size_t                 bin,
                       std::size_t                 order,
                       std::size_t                 channel,
                       std::vector<double>        &node_values,
                       std::vector<double>        &subgrid_array,
                       std::vector<std::size_t>   &shape);

} // namespace pineapfel
