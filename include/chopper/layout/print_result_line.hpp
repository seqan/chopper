#pragma once

#include <fstream>

#include <chopper/data_store.hpp>
#include <chopper/layout/layout.hpp>
#include <chopper/layout/previous_level.hpp>

namespace chopper::layout
{

inline void print_result_line_to(layout::layout::user_bin user_bin, std::ostream & stream)
{
    stream << user_bin.id << '\t';
    for (auto bin : user_bin.previous_TB_indices)
        stream << bin << ';';
    stream << user_bin.storage_TB_id << '\t';
    for (size_t i = 0; i < user_bin.previous_TB_indices.size(); ++i) // number of bins per merged level is 1
        stream << '1' << ';';
    stream << user_bin.number_of_technical_bins;
    stream << '\n';
}

} // namespace chopper::layout
