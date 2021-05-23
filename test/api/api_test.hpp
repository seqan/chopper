#pragma once

#include <seqan3/std/filesystem>
#include <string>

// Generate the full path of a test input file that is provided in the data directory.
static std::filesystem::path data(std::string const & filename)
{
    return std::filesystem::path{std::string{DATADIR}}.concat(filename);
}
