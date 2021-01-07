#pragma once

#include <seqan3/std/charconv>
#include <seqan3/std/ranges>
#include <string>
#include <vector>

#include <seqan3/range/views/to.hpp>

#include <chopper/detail_bin_prefixes.hpp>
#include <chopper/detail_starts_with.hpp>

struct chopper_pack_record
{
    std::string bin_name{};
    int64_t hidx{};
    int64_t lidx{-1};
    std::vector<std::string> filenames{};
    size_t bins{};
    size_t max_size{};

    bool operator==(chopper_pack_record const & other) const
    {
        return std::tie(bin_name, hidx, lidx, filenames, bins, max_size) ==
               std::tie(other.bin_name, other.hidx, other.lidx, other.filenames, other.bins, other.max_size);
    }

    bool operator!=(chopper_pack_record const & other) const
    {
        return std::tie(bin_name, hidx, lidx, filenames, bins, max_size) !=
               std::tie(other.bin_name, other.hidx, other.lidx, other.filenames, other.bins, other.max_size);
    }
};

std::ostream & operator<<(std::ostream & s, chopper_pack_record const & r)
{
    s << "PACK RECORD:"
      << "\n  -> bin_name:" << r.bin_name
      << "\n  -> hidx:" << r.hidx
      << "\n  -> lidx:" << r.lidx
      << "\n  -> filenames: ";
    for (auto const & filename : r.filenames)
        s << filename << ",";
    s << "\n  -> bins:" << r.bins
      << "\n  -> max_size:" << r.max_size
      << '\n';

      return s;
}

auto parse_bin_name(std::string const & name)
{
    int64_t hidx{};
    int64_t lidx{-1};

    if (starts_with(name, merged_bin_prefix))
    {
        // merged bin names look like this: MERGED_BIN_X_Y
        // where X is the high level bin index and Y is the low level bin index
        auto const hend = std::find(name.begin() + merged_bin_prefix_length + 1, name.end(), '_');
        std::string_view::size_type const hlength = (hend - name.begin()) - merged_bin_prefix_length - 1;

        std::string const hidx_str{&name[0] + merged_bin_prefix_length + 1, hlength};
        std::string const lidx_str{hend + 1, name.end()};

        std::from_chars(hidx_str.c_str(), hidx_str.c_str() + hidx_str.size(), hidx);
        std::from_chars(lidx_str.c_str(), lidx_str.c_str() + lidx_str.size(), lidx);
    }
    else
    {
        // split bin names look like this: SPLIT_BIN_X
        // where X is the high level bin index and Y is the low level bin index
        auto const hend = std::find(name.begin() + split_bin_prefix.size() + 1, name.end(), '_');
        std::string_view::size_type const hlength = (hend - name.begin()) - split_bin_prefix.size() - 1;

        std::string const hidx_str{&name[0] + split_bin_prefix.size() + 1, hlength};

        std::from_chars(hidx_str.c_str(), hidx_str.c_str() + hidx_str.size(), hidx);

        lidx = -1; // no low level bin index
    }

    return std::make_tuple(hidx, lidx);
}

auto parse_chopper_pack_line(std::string const & current_line)
{
    chopper_pack_record result{};

    // start parsing
    char const * buffer = current_line.c_str();
    auto field_start = &buffer[0];
    auto field_end = &buffer[0];
    auto const buffer_end = field_start + current_line.size();

    while (field_end != buffer_end && *field_end != '\t') ++field_end;

    result.bin_name = std::string(field_start, field_end);

    std::tie(result.hidx, result.lidx) = parse_bin_name(result.bin_name);

    ++field_end; // skip tab
    field_start = field_end;
    while (field_end != buffer_end && *field_end != '\t') ++field_end;

    std::string filenames_str = std::string(field_start, field_end);
    for (auto && filename : filenames_str | std::views::split(';'))
        result.filenames.push_back((filename | seqan3::views::to<std::string>));

    // read number of technical bins assigned to these files
    ++field_end; // skip tab
    auto res = std::from_chars(field_end, buffer_end, result.bins);

    // read estimated maximum technical bin size
    field_end = res.ptr + 1;
    res = std::from_chars(field_end, buffer_end, result.max_size);

    return result;
}
