#pragma once

#include <sharg/detail/to_string.hpp>
#include <sharg/exceptions.hpp>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <chopper/configuration.hpp>
#include <chopper/sketch/check_filenames.hpp>
#include <chopper/sketch/compute_sketches.hpp>
#include <chopper/sketch/output.hpp>
#include <chopper/sketch/read_data_file.hpp>

namespace chopper::sketch
{

struct dna4_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using sequence_file_type = seqan3::sequence_file_input<dna4_traits,
                                                       seqan3::fields<seqan3::field::seq>,
                                                       seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

inline int
execute(configuration & config, std::vector<std::string> const & filenames, std::vector<sketch::hyperloglog> & sketches)
{
    if (filenames.empty())
        throw sharg::parser_error{
            sharg::detail::to_string("The file ", config.data_file.string(), " appears to be empty.")};

    chopper::sketch::check_filenames(filenames, config);

    if (config.precomputed_files)
    {
        auto hash_minimizer_file = [](auto const & filename)
        {
            std::vector<uint64_t> result{};

            uint64_t hash{};
            char * const hash_data{reinterpret_cast<char *>(&hash)};
            size_t const hash_bytes{sizeof(hash)};

            std::ifstream infile{filename, std::ios::binary};

            while (infile.read(hash_data, hash_bytes))
                result.push_back(hash);

            return result;
        };
        auto hashes = filenames | std::views::transform(hash_minimizer_file);
        sketch::compute_sketches_from_hashes(hashes, config.sketch_bits, config.threads, sketches);
    }
    else
    {
        auto hash_file = [&](auto const & f)
        {
            std::vector<uint64_t> result{};
            sequence_file_type fin{f};

            for (auto && [seq] : fin)
            {
                for (auto hash_value : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.k}))
                    result.push_back(hash_value);
            }

            return result;
        };
        auto hashes = filenames | std::views::transform(hash_file);
        sketch::compute_sketches_from_hashes(hashes, config.sketch_bits, config.threads, sketches);
    }

    if (!config.disable_sketch_output)
    {
        if (!std::filesystem::exists(config.sketch_directory))
            std::filesystem::create_directory(config.sketch_directory);

        assert(filenames.size() == sketches.size());
        for (size_t i = 0; i < filenames.size(); ++i)
            write_sketch_file(filenames[i], sketches[i], config);
    }

    return 0;
}

} // namespace chopper::sketch
