// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <cassert>
#include <filesystem>

#include <cereal/archives/json.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <chopper/prefixes.hpp>

#include <hibf/cereal/path.hpp> // IWYU pragma: keep
#include <hibf/config.hpp>

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
    std::filesystem::path output_filename{"layout.txt"};

    //!\brief Whether the program should determine the best number of IBF bins by doing multiple binning runs.
    bool determine_best_tmax{false};

    //!\brief Whether the programm should compute all binnings up to the given t_max.
    bool force_all_binnings{false};

    //!\brief Whether to print verbose output when computing the statistics when computing the layout.
    bool output_verbose_statistics{false};
    //!\}

    //!\brief The HIBF config which will be used to compute the layout within the HIBF lib.
    seqan::hibf::config hibf_config;

    void read_from(std::istream & stream)
    {
        std::string line;
        std::stringstream config_str;

        while (std::getline(stream, line) && line != chopper::prefix::meta_chopper_config_start)
            ;

        assert(line == chopper::prefix::meta_chopper_config_start);

        while (std::getline(stream, line) && line != chopper::prefix::meta_chopper_config_end)
        {
            assert(line.size() >= 2);
            assert(std::string_view{line}.substr(0, 1) == seqan::hibf::prefix::meta_header);
            config_str << line.substr(1); // remove seqan::hibf::prefix::meta_header
        }

        assert(line == chopper::prefix::meta_chopper_config_end);

        cereal::JSONInputArchive iarchive(config_str);
        iarchive(*this);

        hibf_config.read_from(stream);
    }

    void write_to(std::ostream & stream) const
    {
        // write json file to temprorary string stream with cereal
        std::stringstream config_stream{};
        cereal::JSONOutputArchive output(config_stream); // stream to cout
        output(cereal::make_nvp("chopper_config", *this));

        // write config
        stream << chopper::prefix::meta_chopper_config_start << '\n';
        std::string line;
        while (std::getline(config_stream, line, '\n'))
            stream << seqan::hibf::prefix::meta_header << line << '\n';
        stream << seqan::hibf::prefix::meta_header << "}\n" // last closing bracket isn't written by loop above
               << chopper::prefix::meta_chopper_config_end << '\n';

        hibf_config.write_to(stream);
    }

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
    }
};

} // namespace chopper
