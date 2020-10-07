#pragma once

#include <set>
#include <string>
#include <vector>

#include <seqan3/std/charconv>

#include <chopper/build/build_config.hpp>
#include <chopper/build/data_file_record.hpp>
#include <chopper/detail_starts_with.hpp>

auto parse_line(std::string const & current_line)
{
    // results
    std::string bin_name;
    std::vector<std::string> filenames;
    size_t num_technical_bins;

    // start parsing
    char const * buffer = current_line.c_str();
    auto field_start = &buffer[0];
    auto field_end = &buffer[0];
    auto const buffer_end = field_start + current_line.size();

    while (field_end != buffer_end && *field_end != '\t') ++field_end;

    bin_name = std::string(field_start, field_end);

    ++field_end; // skip tab
    field_start = field_end;
    while (field_end != buffer_end && *field_end != '\t') ++field_end;

    std::string filenames_str = std::string(field_start, field_end);

    // read number of technical bins assigned to these files
    ++field_end; // skip tab
    auto res = std::from_chars(field_end, buffer_end, num_technical_bins);

    for (auto && filename : filenames_str | std::views::split(';'))
        filenames.push_back((filename | seqan3::views::to<std::string>));

    return std::make_tuple(bin_name, filenames, num_technical_bins);
}

auto read_data_file_and_set_high_level_bins(build_config & config)
{
    std::vector<data_file_record> records{};

    std::set<std::string> low_level_name_set{};
    size_t count{};

    std::ifstream binning_file{config.binning_filename};

    if (!binning_file.good() || !binning_file.is_open())
        throw std::logic_error{"Could not open file for reading"};

    std::string current_line;
    std::getline(binning_file, current_line); // skip header line

    while (std::getline(binning_file, current_line))
    {
        auto && [bin_name, filenames, num_technical_bins] = parse_line(current_line);
        records.emplace_back(bin_name, filenames, num_technical_bins);

        if (starts_with(bin_name, merged_bin_prefix))
        {
            std::string prefix(bin_name.begin(),
                               std::find(bin_name.begin() + merged_bin_prefix.size(), bin_name.end(), '_'));
            low_level_name_set.insert(prefix);
        }
        else
        {
            count += num_technical_bins;
        }
    }

    config.high_level_ibf_num_technical_bins = count + low_level_name_set.size();
    assert(config.high_level_ibf_num_technical_bins != 0);

    return records;
};
