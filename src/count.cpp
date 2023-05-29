#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_set>
#include <vector>

#include <sharg/parser.hpp>

#include <chopper/sketch/estimate_kmer_counts.hpp>
#include <chopper/sketch/execute.hpp>
#include <chopper/sketch/read_data_file.hpp>

struct cli_args
{
    std::filesystem::path input_path{};
    std::filesystem::path counts_path{};
    //!\brief The kmer size to hash the input sequences before computing a HyperLogLog sketch from them.
    uint8_t kmer_size{19};
    //!\brief The number of bits the HyperLogLog sketch should use to distribute the values into bins.
    uint8_t sketch_bits{12};
    //!\brief The name for the output directory when writing sketches to disk.
    std::filesystem::path sketch_directory{};
};

int main(int argc, char const * argv[])
{
    sharg::parser parser{"measure_hyperloglog", argc, argv, sharg::update_notifications::off};

    cli_args args{};

    parser.add_option(args.input_path,
                      sharg::config{.short_id = 'i',
                                    .long_id = "input-file",
                                    .description = "Fasta formatted file with sequences.",
                                    .required = true});
    parser.add_option(args.counts_path,
                      sharg::config{.short_id = 'o',
                                    .long_id = "output-file",
                                    .description = "File where the output is written to in tsv format.",
                                    .required = true});
    parser.add_option(args.kmer_size,
                      sharg::config{.short_id = 'k',
                                    .long_id = "kmer-size",
                                    .description = "The size of the k-mers of which the hash values are computed."});

    parser.add_option(
        config.sketch_directory,
        sharg::config{
            .long_id = "output-sketches-to",
            .description =
                "If you supply a directory path with this option, the hyperloglog sketches of your input will be "
                "stored in the respective path; one .hll file per input file.",
            .default_message = "None"});

    parser.add_section("Advanced Options:");
    parser.add_option(
        args.sketch_bits,
        sharg::config{.short_id = 'b',
                      .long_id = "hll-bits",
                      .description = "Adds an integer value which is tested as HyperLogLog bits parameter.",
                      .advances = true});

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[COMMAND LINE INPUT ERROR] " << ext.what() << std::endl;
        return -1;
    }

    std::ofstream fout{args.counts_path};

    std::vector<std::string> filenames{};
    std::vector<size_t> kmer_counts{};
    std::vector<chopper::sketch::hyperloglog> sketches{};

    chopper::sketch::read_data_file(args.input_path, filenames);

    chopper::sketch::execute(config, filenames, sketches);
    chopper::sketch::estimate_kmer_counts(sketches, kmer_counts);

    if (!args.sketch_directory.empty() && !std::filesystem::exists(config.sketch_directory))
        std::filesystem::create_directory(config.sketch_directory);

    for (size_t i = 0; i < filenames.size(); ++i)
    {
        fout << filenames[i] << '\t' << kmer_counts[i] << '\n';
        if (!args.sketch_directory.empty())
            write_sketch_file(filenames[i], sketches[i], args.sketch_directory);
    }
}
