#pragma once

#include <fstream>
#include <string>
#include <vector>

#include <seqan3/utility/views/join_with.hpp>

#include <chopper/configuration.hpp>
#include <chopper/sketch/hyperloglog.hpp>

namespace chopper::sketch
{

inline void write_count_file_line(std::pair<std::string, std::vector<std::string>> const & cluster,
                                  uint64_t const weight,
                                  std::ofstream & fout)
{
    auto & [key, filepaths] = cluster;

    for (auto && arr : filepaths | seqan3::views::join_with(';'))
        fout << arr;

    fout << '\t' << weight << '\t' << key << '\n';
}

inline void write_sketch_file(std::string const & filename,
                              chopper::sketch::hyperloglog const & sketch,
                              std::filesystem::path output_directory)
{
    // For one file in the cluster, the file stem is used with the .hll ending
    std::filesystem::path path = output_directory / std::filesystem::path(filename).stem();
    path += ".hll";
    std::ofstream hll_fout(path, std::ios::binary);
    sketch.dump(hll_fout);
}

} // namespace chopper::sketch
