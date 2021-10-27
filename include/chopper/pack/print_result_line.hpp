#pragma once

#include <fstream>

#include <chopper/pack/pack_data.hpp>
#include <chopper/pack/previous_level.hpp>

namespace chopper::pack
{

inline void print_result_line(pack_data const & data,
                              size_t const index,
                              size_t const bin_id,
                              size_t const number_of_bins)
{
    bool const is_top_level = data.previous.empty();

    *data.output_buffer << data.filenames[index]
                        << '\t' << data.previous.bin_indices << (is_top_level ? "" : ";") << bin_id
                        << '\t' << data.previous.num_of_bins << (is_top_level ? "" : ";") << number_of_bins
                        << '\n';
}

inline void print_debug_line(pack_data const & data,
                             size_t const index,
                             size_t const bin_id,
                             size_t const number_of_bins,
                             size_t const average_bin_size,
                             size_t const optimal_score,
                             size_t const num_technical_bins)
{
    bool const is_top_level = data.previous.empty();

    assert(number_of_bins > 0);
    double const correction = data.fp_correction[number_of_bins];

    *data.output_buffer << data.filenames[index]
                        << '\t' << data.previous.bin_indices << (is_top_level ? "" : ";") << bin_id
                        << '\t' << data.previous.num_of_bins << (is_top_level ? "" : ";") << number_of_bins
                        << '\t' << data.previous.estimated_sizes << (is_top_level ? "" : ";") << average_bin_size
                        << '\t' << data.previous.optimal_score << (is_top_level ? "" : ";") << optimal_score
                        << '\t' << data.previous.correction << (is_top_level ? "" : ";") << correction
                        << '\t' << data.previous.tmax << (is_top_level ? "" : ";") << num_technical_bins
                        << '\n';
}

} // namespace chopper::pack
