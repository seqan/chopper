#pragma once

#include <seqan3/std/charconv>
#include <seqan3/std/ranges>
#include <string>
#include <vector>

#include <seqan3/utility/views/to.hpp>

#include <chopper/detail_bin_prefixes.hpp>
#include <chopper/detail_starts_with.hpp>

struct chopper_pack_record
{
    std::vector<std::string> filenames{};
    std::vector<size_t> bin_indices{};
    std::vector<size_t> number_of_bins{};
    std::vector<size_t> estimated_sizes{};

    bool operator==(chopper_pack_record const & other) const
    {
        return std::tie(filenames, bin_indices, number_of_bins, estimated_sizes) ==
               std::tie(other.filenames, other.bin_indices, other.number_of_bins, other.estimated_sizes);
    }

    bool operator!=(chopper_pack_record const & other) const
    {
        return std::tie(filenames, bin_indices, number_of_bins, estimated_sizes) !=
               std::tie(other.filenames, other.bin_indices, other.number_of_bins, other.estimated_sizes);
    }
};

inline std::ostream & operator<<(std::ostream & s, chopper_pack_record const & r)
{
    s << "PACK RECORD:"
      << "\n  -> filenames: ";
    for (auto const & filename : r.filenames)
        s << filename << ",";
    s << "\n  -> bin_indices: ";
    for (auto const & index : r.bin_indices)
        s << index << ",";
    s  << "\n  -> number_of_bins: ";
    for (auto const & num : r.number_of_bins)
        s << num << ",";
    s << "\n  -> estimated_sizes: ";
    for (auto const & estimate : r.estimated_sizes)
        s << estimate << ",";
    s << '\n';

      return s;
}

inline auto parse_chopper_pack_line(std::string const & current_line)
{
    chopper_pack_record result{};

    // initialize parsing
    char const * buffer = current_line.c_str();
    auto field_start = &buffer[0];
    auto field_end = &buffer[0];
    auto const buffer_end = field_start + current_line.size();
    while (field_end != buffer_end && *field_end != '\t') ++field_end;

    // parse filenames
    std::string filenames_str = std::string(field_start, field_end);
    for (auto && filename : filenames_str | std::views::split(';'))
        result.filenames.push_back((filename | seqan3::views::to<std::string>));

    size_t tmp; // temporary size_t

    do // read bin_indices
    {
        ++field_end; // skip tab or ;
        auto res = std::from_chars(field_end, buffer_end, tmp);
        field_end = res.ptr;
        result.bin_indices.push_back(tmp);
    } while (field_end != buffer_end && *field_end != '\t');

    do // read number of technical bins
    {
        ++field_end; // skip tab or ;
        auto res = std::from_chars(field_end, buffer_end, tmp);
        field_end = res.ptr;
        result.number_of_bins.push_back(tmp);
    } while (field_end != buffer_end && *field_end != '\t');

    do // read estimated maximum technical bin size
    {
        ++field_end; // skip tab or ;
        auto res = std::from_chars(field_end, buffer_end, tmp);
        field_end = res.ptr;
        result.estimated_sizes.push_back(tmp);
    } while (field_end != buffer_end && *field_end != '\n');

    return result;
}
