#pragma once

#include <filesystem>
#include <fstream>
#include <future>
#include <thread>

#include <robin_hood.h>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/views/async_input_buffer.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/range/to.hpp>

#include <chopper/configuration.hpp>
#include <chopper/count/output.hpp>
#include <chopper/sketch/hyperloglog.hpp>

namespace chopper::count
{

struct mytraits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using sequence_file_type = seqan3::sequence_file_input<mytraits,
                                                       seqan3::fields<seqan3::field::seq>,
                                                       seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

inline void process_sequence_files(std::vector<std::string> const & filenames,
                                   configuration const & config,
                                   sketch::hyperloglog & sketch)
{
    for (auto const & filename : filenames)
        for (auto && [seq] : sequence_file_type{filename})
            for (auto && k_hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
                sketch.add(reinterpret_cast<char*>(&k_hash), sizeof(k_hash));
}

inline void process_minimizer_files(std::vector<std::string> const & filenames, sketch::hyperloglog & sketch)
{
    // temporary variables when .minimizer files are read
    uint64_t hash{};
    char * const hash_data{reinterpret_cast<char*>(&hash)};
    size_t const hash_bytes{sizeof(hash)};

    // read files
    for (auto const & filename : filenames)
    {
        std::ifstream infile{filename, std::ios::binary};

        while (infile.read(hash_data, hash_bytes))
            sketch.add(hash_data, hash_bytes);
    }
}

inline void count_kmers(robin_hood::unordered_map<std::string, std::vector<std::string>> const & filename_clusters,
                        configuration const & config)
{
   // output file
    std::ofstream fout{config.count_filename};

    if (!fout.good())
        throw std::runtime_error{"Could not open file" + config.count_filename.string() + " for reading."};

    // create the hll dir if it doesn't already exist
    if (!config.disable_sketch_output)
        std::filesystem::create_directory(config.sketch_directory);

    // copy filename clusters to vector
    std::vector<std::pair<std::string, std::vector<std::string>>> cluster_vector{};
    for (auto const & cluster : filename_clusters)
        cluster_vector.emplace_back(cluster.first, cluster.second);

    #pragma omp parallel for schedule(static) num_threads(config.threads)
    for (size_t i = 0; i < cluster_vector.size(); ++i)
    {
        chopper::sketch::hyperloglog sketch(config.sketch_bits);

        if (config.precomputed_files)
            process_minimizer_files(cluster_vector[i].second, sketch);
        else
            process_sequence_files(cluster_vector[i].second, config, sketch);

        // print either the exact or the approximate count, depending on exclusively_hlls
        uint64_t const weight = sketch.estimate();

        #pragma omp critical
        write_count_file_line(cluster_vector[i], weight, fout);

        if (!config.disable_sketch_output)
            write_sketch_file(cluster_vector[i], sketch, config);
    }
}

} // namespace chopper::count
