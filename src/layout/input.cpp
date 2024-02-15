// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <cassert>
#include <charconv>
#include <istream>
#include <ranges>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <chopper/configuration.hpp>
#include <chopper/layout/input.hpp>
#include <chopper/prefixes.hpp>

#include <hibf/layout/layout.hpp>

namespace chopper::layout
{

std::vector<std::vector<std::string>> read_filenames_from(std::istream & stream)
{
    std::vector<std::vector<std::string>> filenames{};
    std::string line;

    while (std::getline(stream, line) && line != chopper::prefix::meta_chopper_user_bins_start)
        ;

    assert(line == chopper::prefix::meta_chopper_user_bins_start);

#ifndef NDEBUG
    size_t counter{};
#endif
    while (std::getline(stream, line) && line != chopper::prefix::meta_chopper_user_bins_end)
    {
        assert(line.size() >= 2);
        assert(std::string_view{line}.substr(0, 1) == seqan::hibf::prefix::meta_header);

        // @0 file1.fa file2.fa
        auto const bin_idx_pos = line.find(' ');
        assert(bin_idx_pos != std::string::npos);

#ifndef NDEBUG
        size_t bin_idx{};
        std::from_chars(line.data() + 1, line.data() + bin_idx_pos, bin_idx);
        assert(bin_idx == counter++);
#endif

        filenames.emplace_back();
        std::string_view const filename_str{line.begin() + bin_idx_pos + 1, line.end()};
        for (auto const && filename : std::views::split(filename_str, ' '))
        {
            auto common_view = std::views::common(filename);
            filenames.back().emplace_back(common_view.begin(), common_view.end());
        }
    }

    assert(line == chopper::prefix::meta_chopper_user_bins_end);

    return filenames;
}

std::tuple<std::vector<std::vector<std::string>>, configuration, std::vector<seqan::hibf::layout::layout>>
read_layouts_file(std::istream & stream)
{
    std::vector<std::vector<std::string>> filenames = chopper::layout::read_filenames_from(stream);
    chopper::configuration chopper_config;
    chopper_config.read_from(stream);
    std::vector<seqan::hibf::layout::layout> layouts;
    while (stream.good())
    {
        seqan::hibf::layout::layout hibf_layout{};
        hibf_layout.read_from(stream);
        layouts.push_back(std::move(hibf_layout));
    }
    return std::make_tuple(std::move(filenames), std::move(chopper_config), std::move(layouts));
}

} // namespace chopper::layout
