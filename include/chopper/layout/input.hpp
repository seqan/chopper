// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <cassert>
#include <charconv>
#include <ranges>

#include <cereal/archives/json.hpp>

#include <chopper/configuration.hpp>

#include <hibf/detail/prefixes.hpp>

namespace chopper::layout
{

inline void read_config_from(configuration & config, std::istream & stream)
{
    std::string line;
    std::stringstream config_str;

    while (std::getline(stream, line) && line != "@CHOPPER_CONFIG")
        ;

    assert(line == "@CHOPPER_CONFIG");

    // TODO ##CONFIG: as prefix
    while (std::getline(stream, line) && line != "@CHOPPER_CONFIG_END")
    {
        assert(line.size() >= 2);
        assert(std::string_view{line}.substr(0, 1) == seqan::hibf::prefix::meta_header);
        config_str << line.substr(1); // remove seqan::hibf::prefix::meta_header
    }

    assert(line == "@CHOPPER_CONFIG_END");

    cereal::JSONInputArchive iarchive(config_str);
    iarchive(config);
}

inline std::vector<std::vector<std::string>> read_filenames_from(std::istream & stream)
{
    std::vector<std::vector<std::string>> filenames{};
    std::string line;

    while (std::getline(stream, line) && line != "@CHOPPER_USER_BINS")
        ;

    assert(line == "@CHOPPER_USER_BINS");

#ifndef NDEBUG
    size_t counter{};
#endif
    // TODO ##CONFIG: as prefix
    while (std::getline(stream, line) && line != "@CHOPPER_USER_BINS_END")
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

    assert(line == "@CHOPPER_USER_BINS_END");

    return filenames;
}

} // namespace chopper::layout
