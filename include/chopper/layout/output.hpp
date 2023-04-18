#pragma once

#include <cereal/archives/json.hpp>

#include <chopper/configuration.hpp>
#include <chopper/layout/layout.hpp>
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

inline void write_layout_header_to(layout const & hibf_layout, size_t const max_hibf_id, std::ostream & stream)
{
    stream << prefix::first_header_line << " max_bin_id:" << max_hibf_id << '\n';
    for (auto const & max_bin : hibf_layout.max_bins)
        stream << max_bin << '\n';
}

inline void write_user_bin_line_to(layout::user_bin const & object,
                                   std::vector<std::string> const & filenames,
                                   std::ostream & stream)
{
    stream << filenames[object.idx] << '\t';
    for (auto bin : object.previous_TB_indices)
        stream << bin << ';';
    stream << object.storage_TB_id << '\t';
    for ([[maybe_unused]] auto && elem : object.previous_TB_indices) // number of bins per merged level is 1
        stream << "1;";
    stream << object.number_of_technical_bins;
    stream << '\n';
}

inline void
write_layout_content_to(layout const & hibf_layout, std::vector<std::string> const & filenames, std::ostream & stream)
{
    stream << prefix::header << "FILES\tBIN_INDICES\tNUMBER_OF_BINS\n";
    for (auto const & user_bin : hibf_layout.user_bins)
        write_user_bin_line_to(user_bin, filenames, stream);
}

} // namespace chopper::layout
