#pragma once

#include <seqan3/std/charconv>
#include <seqan3/std/ranges>
#include <string>
#include <vector>

#include <seqan3/utility/views/to.hpp>

#include <chopper/detail_bin_prefixes.hpp>
#include <chopper/detail_starts_with.hpp>
#include <chopper/detail_string_view.hpp>

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

inline chopper_pack_record parse_chopper_pack_line(std::string const & current_line)
{
    chopper_pack_record result{};

    // initialize parsing
    std::string_view const buffer{current_line};
    auto const buffer_end{buffer.end()};
    auto field_end = buffer.begin();
    while (field_end != buffer_end && *field_end != '\t') ++field_end;

    // parse filenames
    std::string_view const filenames{chopper::detail::string_view(buffer.begin(), field_end)};
    for (auto && filename : filenames | std::views::split(';'))
    {
        auto const common_view = filename | std::views::common;
        result.filenames.emplace_back(common_view.begin(), common_view.end());
    }

    size_t tmp{};

    do // read bin_indices
    {
        ++field_end; // skip tab or ;
        field_end = std::from_chars(field_end, buffer_end, tmp).ptr;
        result.bin_indices.push_back(tmp);
    } while (field_end != buffer_end && *field_end != '\t');

    do // read number of technical bins
    {
        ++field_end; // skip tab or ;
        field_end = std::from_chars(field_end, buffer_end, tmp).ptr;
        result.number_of_bins.push_back(tmp);
    } while (field_end != buffer_end && *field_end != '\t');

    return result;
}
