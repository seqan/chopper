#pragma once

#include <seqan3/std/filesystem>
#include <fstream>
#include <string>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/tmp_filename.hpp>

// Generate the full path of a test input file that is provided in the data directory.
[[maybe_unused]] static std::filesystem::path data(std::string const & filename)
{
    return std::filesystem::path{std::string{DATADIR}}.concat(filename);
}

// Returns contents of a file as string.
[[maybe_unused]] static std::string const string_from_file(std::filesystem::path const & path,
                                                           std::ios_base::openmode const mode = std::ios_base::in)
{
    std::ifstream file_stream(path, mode);
    if (!file_stream.is_open())
        throw std::logic_error{"Cannot open " + path.string()};
    std::stringstream file_buffer;
    file_buffer << file_stream.rdbuf();
    return {file_buffer.str()};
}
