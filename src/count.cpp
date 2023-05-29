#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <sharg/parser.hpp>

#include <chopper/configuration.hpp>
#include <chopper/sketch/estimate_kmer_counts.hpp>
#include <chopper/sketch/execute.hpp>
#include <chopper/sketch/read_data_file.hpp>

int main(int argc, char const * argv[])
{
    sharg::parser parser{"chopper_count", argc, argv, sharg::update_notifications::off};

    chopper::configuration config;
    std::string counts_path;

    parser.add_option(config.data_file,
                      sharg::config{.short_id = 'i',
                                    .long_id = "input-file",
                                    .description = "Fasta formatted file with sequences.",
                                    .required = true});
    parser.add_option(counts_path,
                      sharg::config{.short_id = 'o',
                                    .long_id = "output-file",
                                    .description = "File where the output is written to in tsv format.",
                                    .required = true});
    parser.add_option(config.k,
                      sharg::config{.short_id = 'k',
                                    .long_id = "kmer-size",
                                    .description = "The size of the k-mers of which the hash values are computed."});
    parser.add_option(config.threads,
                      sharg::config{.short_id = 't',
                                    .long_id = "threads",
                                    .description = "The number of threads to use for sketching."});
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
        config.sketch_bits,
        sharg::config{.short_id = 'b',
                      .long_id = "hll-bits",
                      .description = "Adds an integer value which is tested as HyperLogLog bits parameter.",
                      .advanced = true});

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[COMMAND LINE INPUT ERROR] " << ext.what() << std::endl;
        return -1;
    }

    config.disable_sketch_output = !parser.is_option_set("output-sketches-to");

    std::ofstream fout{counts_path};

    std::vector<std::string> filenames{};
    std::vector<size_t> kmer_counts{};
    std::vector<chopper::sketch::hyperloglog> sketches{};

    chopper::sketch::read_data_file(config.data_file, filenames);

    chopper::sketch::execute(config, filenames, sketches);
    chopper::sketch::estimate_kmer_counts(sketches, kmer_counts);

    for (size_t i = 0; i < filenames.size(); ++i)
        fout << filenames[i] << '\t' << kmer_counts[i] << '\n';
}
