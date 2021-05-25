#pragma once

#include <seqan3/std/filesystem>
#include <string>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/tmp_filename.hpp>

// Generate the full path of a test input file that is provided in the data directory.
static std::filesystem::path data(std::string const & filename)
{
    return std::filesystem::path{std::string{DATADIR}}.concat(filename);
}
