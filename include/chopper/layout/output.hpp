// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <cereal/archives/json.hpp>

#include <chopper/configuration.hpp>
#include <chopper/prefixes.hpp>

#include <hibf/detail/layout/layout.hpp>

namespace chopper::layout
{

inline void write_config_to(configuration const & config, std::ostream & stream)
{
    // write json file to temprorary string stream with cereal
    std::stringstream config_stream{};
    cereal::JSONOutputArchive output(config_stream); // stream to cout
    output(cereal::make_nvp("chopper_config", config));

    // write config
    stream << seqan::hibf::prefix::meta_header << "CHOPPER_CONFIG\n";
    std::string line;
    while (std::getline(config_stream, line, '\n'))
        stream << seqan::hibf::prefix::meta_header << line << '\n';
    stream << seqan::hibf::prefix::meta_header << "}\n" // last closing bracket isn't written by loop above
           << seqan::hibf::prefix::meta_header << "CHOPPER_CONFIG_END\n";

    config.hibf_config.write_to(stream);
}

inline void write_user_bins_to(std::vector<std::string> const & filenames, std::ostream & stream)
{
    stream << seqan::hibf::prefix::meta_header << "CHOPPER_USER_BINS\n";
    size_t counter{};
    for (auto const & filename : filenames)
        stream << seqan::hibf::prefix::meta_header << counter++ << ' ' << filename << '\n';
    stream << seqan::hibf::prefix::meta_header << "CHOPPER_USER_BINS_END\n";
}

} // namespace chopper::layout
