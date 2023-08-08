// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <filesystem>

#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <hibf/config.hpp>
#include <hibf/detail/cereal/path.hpp> // IWYU pragma: keep

namespace chopper
{

struct configuration
{
    /*!\name General Configuration
     * \{
     */
    //!\brief The input file to chopper. Should contain one file path per line.
    std::filesystem::path data_file;

    //!\brief Internal parameter that triggers some verbose debug output.
    bool debug{false};
    //!\}

    /*!\name Configuration of size estimates (chopper::count)
     * \{
     */
    //!\brief The name for the output directory when writing sketches to disk.
    std::filesystem::path sketch_directory{};

    //!\brief The kmer size to hash the input sequences before computing a HyperLogLog sketch from them.
    uint8_t k{19};

    //!\brief Do not write the sketches into a dedicated directory.
    bool disable_sketch_output{false};

    //!\brief Whether the input files are precomputed files (.minimiser) instead of sequence files.
    bool precomputed_files{false};
    //!\}

    /*!\name General Configuration
     * \{
     */
    //!\brief The name of the layout file to write.
    std::filesystem::path output_filename{"binning.out"};

    //!\brief Whether the program should determine the best number of IBF bins by doing multiple binning runs.
    bool determine_best_tmax{false};

    //!\brief Whether the programm should compute all binnings up to the given t_max.
    bool force_all_binnings{false};

    //!\brief Whether to print verbose output when computing the statistics when computing the layout.
    bool output_verbose_statistics{false};
    //!\}

    //!\brief The HIBF config which will be used to compute the layout within the HIBF lib.
    hibf::config hibf_config;

private:
    friend class cereal::access;

    template <typename archive_t>
    void serialize(archive_t & archive)
    {
        uint32_t version{2};
        archive(CEREAL_NVP(version));

        archive(CEREAL_NVP(data_file));
        archive(CEREAL_NVP(debug));
        archive(CEREAL_NVP(sketch_directory));
        archive(CEREAL_NVP(k));
        archive(CEREAL_NVP(disable_sketch_output));
        archive(CEREAL_NVP(precomputed_files));

        archive(CEREAL_NVP(output_filename));
        archive(CEREAL_NVP(determine_best_tmax));
        archive(CEREAL_NVP(force_all_binnings));

        archive(CEREAL_NVP(hibf_config));
    }
};

} // namespace chopper
