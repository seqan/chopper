#pragma once

#include <string>
#include <vector>

#include <seqan3/std/charconv>
#include <seqan3/std/ranges>

#include <seqan3/range/views/to.hpp>

struct data_file_record
{
    std::string bin_name{};
    std::vector<std::string> filenames{};
    size_t bins{};
};

auto parse_binning_line(std::string const & current_line)
{
    data_file_record result{};

    // start parsing
    char const * buffer = current_line.c_str();
    auto field_start = &buffer[0];
    auto field_end = &buffer[0];
    auto const buffer_end = field_start + current_line.size();

    while (field_end != buffer_end && *field_end != '\t') ++field_end;

    result.bin_name = std::string(field_start, field_end);

    ++field_end; // skip tab
    field_start = field_end;
    while (field_end != buffer_end && *field_end != '\t') ++field_end;

    std::string filenames_str = std::string(field_start, field_end);
    for (auto && filename : filenames_str | std::views::split(';'))
        result.filenames.push_back((filename | seqan3::views::to<std::string>));

    // read number of technical bins assigned to these files
    ++field_end; // skip tab
    auto res = std::from_chars(field_end, buffer_end, result.bins);

    return result;
}
