#pragma once

#include <string>
#include <tuple>

#include <seqan3/std/charconv>

#include <chopper/build/region.hpp>

auto parse_chopper_split_file_line(std::string const & line)
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

    region reg{};

    // read begin
    assert(*field_end == '\t');
    ++field_end;
    auto res = std::from_chars(field_end, buffer_end, reg.begin);
    field_end = res.ptr;

    // read end
    assert(*field_end == '\t');
    ++field_end;
    res = std::from_chars(field_end, buffer_end, reg.end);
    field_end = res.ptr;

    // read hibf bin index
    assert(*field_end == '\t');
    ++field_end;
    res = std::from_chars(field_end, buffer_end, reg.hidx);
    field_end = res.ptr;

    // read libf bin index
    assert(*field_end == '\t');
    ++field_end;

    if (*field_end == '-')
        reg.lidx = -1;
    else
        res = std::from_chars(field_end, buffer_end, reg.lidx);

    return std::make_tuple(std::move(filename), std::move(id), std::move(reg));
}
