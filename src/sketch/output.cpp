// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <chopper/sketch/output.hpp>

#include <hibf/contrib/std/join_with_view.hpp>

namespace chopper::sketch
{

void write_count_file_line(std::pair<std::string, std::vector<std::string>> const & cluster,
                           uint64_t const weight,
                           std::ofstream & fout)
{
    auto & [key, filepaths] = cluster;

    for (auto && arr : filepaths | seqan::stl::views::join_with(';'))
        fout << arr;

    fout << '\t' << weight << '\t' << key << '\n';
}

void write_sketch_file(std::string const & filename,
                       seqan::hibf::sketch::hyperloglog const & sketch,
                       configuration const & config)
{
    // For one file in the cluster, the file stem is used with the .hll ending
    std::filesystem::path path = config.sketch_directory / std::filesystem::path(filename).stem();
    path += ".hll";
    std::ofstream hll_fout(path, std::ios::binary);
    sketch.store(hll_fout);
}

} // namespace chopper::sketch
