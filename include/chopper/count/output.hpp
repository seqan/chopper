#pragma once

#include <fstream>
#include <string>
#include <vector>

#include <seqan3/utility/views/join_with.hpp>

#include <chopper/configuration.hpp>
#include <chopper/sketch/hyperloglog.hpp>

namespace chopper::count
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

inline void write_sketch_file(std::pair<std::string, std::vector<std::string>> const & cluster,
                              chopper::sketch::hyperloglog const & sketch,
                              configuration const & config)
{
    auto & [key, filepaths] = cluster;
    // For more than one file in the cluster, Felix doesn't know how to name the file
    // and what exactly is supposed to happen.
    if (filepaths.size() != 1)
        throw std::runtime_error("This mode is not implemented yet for multiple files grouped together.");

    // For one file in the cluster, the file stem is used with the .hll ending
    std::filesystem::path path = config.sketch_directory / std::filesystem::path(filepaths[0]).stem();
    path += ".hll";
    std::ofstream hll_fout(path, std::ios::binary);
    sketch.dump(hll_fout);
}


} // namespace chopper::count
