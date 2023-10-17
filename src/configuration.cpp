// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <chopper/configuration.hpp>

namespace chopper
{

void configuration::read_from(std::istream & stream)
{
    std::string line;
    std::stringstream config_str;

    while (std::getline(stream, line) && line != chopper::prefix::meta_chopper_config_start)
        ;

    assert(line == chopper::prefix::meta_chopper_config_start);

    while (std::getline(stream, line) && line != chopper::prefix::meta_chopper_config_end)
    {
        assert(line.size() >= 2);
        assert(std::string_view{line}.substr(0, 1) == seqan::hibf::prefix::meta_header);
        config_str << line.substr(1); // remove seqan::hibf::prefix::meta_header
    }

    assert(line == chopper::prefix::meta_chopper_config_end);

    cereal::JSONInputArchive iarchive(config_str);
    iarchive(*this);

    hibf_config.read_from(stream);
}

void configuration::write_to(std::ostream & stream) const
{
    // write json file to temprorary string stream with cereal
    std::stringstream config_stream{};
    cereal::JSONOutputArchive output(config_stream); // stream to cout
    output(cereal::make_nvp("chopper_config", *this));

    // write config
    stream << chopper::prefix::meta_chopper_config_start << '\n';
    std::string line;
    while (std::getline(config_stream, line, '\n'))
        stream << seqan::hibf::prefix::meta_header << line << '\n';
    stream << seqan::hibf::prefix::meta_header << "}\n" // last closing bracket isn't written by loop above
           << chopper::prefix::meta_chopper_config_end << '\n';

    hibf_config.write_to(stream);
}

} // namespace chopper
