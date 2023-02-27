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
#include <chopper/data_store.hpp>
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

inline void
process_sequence_file(std::string const & filename, configuration const & config, sketch::hyperloglog & sketch)
{
    for (auto && [seq] : sequence_file_type{filename})
        for (auto && k_hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
            sketch.add(reinterpret_cast<char *>(&k_hash), sizeof(k_hash));
}

inline void process_minimiser_file(std::string const & filename, sketch::hyperloglog & sketch)
{
    // temporary variables when .minimiser files are read
    uint64_t hash{};
    char * const hash_data{reinterpret_cast<char *>(&hash)};
    size_t const hash_bytes{sizeof(hash)};

    // read files
    std::ifstream infile{filename, std::ios::binary};

    while (infile.read(hash_data, hash_bytes))
        sketch.add(hash_data, hash_bytes);
}

inline void count_kmers(configuration const & config, data_store & data)
{
    // create the hll dir if it doesn't already exist
    if (!config.disable_sketch_output)
        std::filesystem::create_directory(config.sketch_directory);

    data.sketches.resize(data.filenames.size());

#pragma omp parallel for schedule(static) num_threads(config.threads)
    for (size_t i = 0; i < data.filenames.size(); ++i)
    {
        chopper::sketch::hyperloglog sketch(config.sketch_bits);

        if (config.precomputed_files)
            process_minimiser_file(data.filenames[i], sketch);
        else
            process_sequence_file(data.filenames[i], config, sketch);

#pragma omp critical
        data.sketches[i] = sketch;

        if (!config.disable_sketch_output)
            write_sketch_file(data.filenames[i], sketch, config);
    }
}

} // namespace chopper::count
