#pragma once

#include <string>
#include <tuple>

#include <seqan3/std/charconv>

auto parse_traversal_file_line(std::string const & line)
{
    char const * buffer = line.c_str();
    auto field_start = &buffer[0];
    auto field_end = &buffer[0];
    auto const buffer_end = field_start + line.size();

    while (field_end != buffer_end && *field_end != '\t') ++field_end;

    std::string const filename(field_start, field_end);

    assert(*field_end == '\t');
    ++field_end; // skip tab
    field_start = field_end;
    while (field_end != buffer_end && *field_end != '\t') ++field_end;

    std::string const id(field_start, field_end);

    // read begin
    size_t begin{};
    assert(*field_end == '\t');
    ++field_end;
    auto res = std::from_chars(field_end, buffer_end, begin);
    field_end = res.ptr;

    // read end
    size_t end{};
    assert(*field_end == '\t');
    ++field_end;
    res = std::from_chars(field_end, buffer_end, end);
    field_end = res.ptr;

    // read bin index
    size_t bin_index{};
    assert(*field_end == '\t');
    ++field_end;
    res = std::from_chars(field_end, buffer_end, bin_index);

    return std::make_tuple(std::move(filename), std::move(id), begin, end, bin_index);
}
