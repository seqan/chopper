#pragma once

#include <cereal/archives/json.hpp>

#include <chopper/configuration.hpp>
#include <chopper/prefixes.hpp>

namespace chopper::layout
{

inline void write_config_to(configuration const & config, std::ostream & stream)
{
    // write json file to temprorary string stream with cereal
    std::stringstream config_stream{};
    cereal::JSONOutputArchive output(config_stream); // stream to cout
    output(cereal::make_nvp("config", config));

    // write config
    stream << prefix::header << prefix::header_config << "CONFIG:\n";
    std::string line;
    while (std::getline(config_stream, line, '\n'))
        stream << prefix::header << prefix::header_config << line << '\n';
    stream << prefix::header << prefix::header_config << "}\n" // last closing bracket isn't written by loop above
           << prefix::header << prefix::header_config << "ENDCONFIG\n";
}

inline void write_layout_header_to(configuration const & config, size_t const max_hibf_id, std::string_view const header, std::ostream & stream)
{
    write_config_to(config, stream);
    stream << prefix::header << prefix::high_level << " max_bin_id:" << max_hibf_id << '\n';
    stream << header;
}

} // namespace chopper::layout
