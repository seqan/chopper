#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <chopper/count/execute.hpp>
#include <chopper/detail_apply_prefix.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/sketch/estimate_kmer_counts.hpp>

#include "../api_test.hpp"

TEST(execute_estimation_test, few_ubs)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};
    std::filesystem::path const stats_file{layout_file.string() + ".stats"};

    chopper::configuration config{};
    config.tmax = 64;
    config.determine_best_tmax = true;
    config.disable_sketch_output = true;
    config.output_filename = layout_file;
    config.disable_estimate_union = true; // also disables rearrangement

    std::stringstream output_buffer;
    std::stringstream header_buffer;

    chopper::data_store store{.false_positive_rate = config.false_positive_rate,
                              .output_buffer = &output_buffer,
                              .header_buffer = &header_buffer,
                              .filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7"},
                              .kmer_counts = {500, 1000, 500, 500, 500, 500, 500, 500}};

    chopper::layout::execute(config, store);

    ASSERT_TRUE(std::filesystem::exists(stats_file));

    std::string const written_file{string_from_file(stats_file)};

    EXPECT_EQ(written_file,
              R"expected_cout(## ### Parameters ###
## number of user bins = 8
## number of hash functions = 2
## false positive rate = 0.05
## ### Notation ###
## X-IBF = An IBF with X number of bins.
## X-HIBF = An HIBF with tmax = X, e.g a maximum of X technical bins on each level.
## ### Column Description ###
## tmax : The maximum number of technical bin on each level
## c_tmax : The technical extra cost of querying an tmax-IBF, compared to 64-IBF
## l_tmax : The estimated query cost for an tmax-HIBF, compared to an 64-HIBF
## m_tmax : The estimated memory consumption for an tmax-HIBF, compared to an 64-HIBF
## (l*m)_tmax : Computed by l_tmax * m_tmax
## size : The expected total size of an tmax-HIBF
# tmax	c_tmax	l_tmax	m_tmax	(l*m)_tmax	size
64	1.00	1.00	1.00	1.00	14.6KiB
# Best t_max (regarding expected query runtime): 64
)expected_cout");
}

TEST(execute_estimation_test, many_ubs)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};
    std::filesystem::path const stats_file{layout_file.string() + ".stats"};

    std::vector<std::string> many_filenames;
    std::vector<size_t> many_kmer_counts;

    // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
    for (size_t i{0}; i < 96u; ++i)
    {
        many_filenames.push_back(seqan3::detail::to_string("seq", i));
        many_kmer_counts.push_back(100 * ((i + 20) / 20));
    }

    chopper::configuration config{};
    config.tmax = 1024;
    config.determine_best_tmax = true;
    config.output_filename = layout_file;
    config.disable_sketch_output = true;
    config.disable_estimate_union = true; // also disables rearrangement

    std::stringstream output_buffer;
    std::stringstream header_buffer;

    chopper::data_store data{.false_positive_rate = config.false_positive_rate,
                             .output_buffer = &output_buffer,
                             .header_buffer = &header_buffer,
                             .filenames = many_filenames,
                             .kmer_counts = many_kmer_counts};

    chopper::layout::execute(config, data);

    ASSERT_TRUE(std::filesystem::exists(stats_file));

    std::string const written_file{string_from_file(stats_file)};

    EXPECT_EQ(written_file,
              R"expected_cout(## ### Parameters ###
## number of user bins = 96
## number of hash functions = 2
## false positive rate = 0.05
## ### Notation ###
## X-IBF = An IBF with X number of bins.
## X-HIBF = An HIBF with tmax = X, e.g a maximum of X technical bins on each level.
## ### Column Description ###
## tmax : The maximum number of technical bin on each level
## c_tmax : The technical extra cost of querying an tmax-IBF, compared to 64-IBF
## l_tmax : The estimated query cost for an tmax-HIBF, compared to an 64-HIBF
## m_tmax : The estimated memory consumption for an tmax-HIBF, compared to an 64-HIBF
## (l*m)_tmax : Computed by l_tmax * m_tmax
## size : The expected total size of an tmax-HIBF
# tmax	c_tmax	l_tmax	m_tmax	(l*m)_tmax	size
64	1.00	1.26	1.00	1.26	73.3KiB
128	1.22	1.25	0.66	0.82	48.4KiB
256	1.33	1.33	0.74	0.99	54.3KiB
# Best t_max (regarding expected query runtime): 128
)expected_cout");

    std::string const layout_string{string_from_file(layout_file)};
    EXPECT_NE(layout_string.find("\"tmax\": 128,"), std::string::npos);
}

TEST(execute_estimation_test, many_ubs_force_all)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};
    std::filesystem::path const stats_file{layout_file.string() + ".stats"};

    std::vector<std::string> many_filenames;
    std::vector<size_t> many_kmer_counts;

    // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
    for (size_t i{0}; i < 96u; ++i)
    {
        many_filenames.push_back(seqan3::detail::to_string("seq", i));
        many_kmer_counts.push_back(100 * ((i + 20) / 20));
    }

    chopper::configuration config{};
    config.tmax = 256;
    config.determine_best_tmax = true;
    config.force_all_binnings = true;
    config.disable_sketch_output = true;
    config.output_filename = layout_file;
    config.disable_estimate_union = true; // also disables rearrangement

    std::stringstream output_buffer;
    std::stringstream header_buffer;

    chopper::data_store data{.false_positive_rate = config.false_positive_rate,
                             .output_buffer = &output_buffer,
                             .header_buffer = &header_buffer,
                             .filenames = many_filenames,
                             .kmer_counts = many_kmer_counts};

    chopper::layout::execute(config, data);

    ASSERT_TRUE(std::filesystem::exists(stats_file));

    std::string const written_file{string_from_file(stats_file)};

    EXPECT_EQ(written_file,
              R"expected_cout(## ### Parameters ###
## number of user bins = 96
## number of hash functions = 2
## false positive rate = 0.05
## ### Notation ###
## X-IBF = An IBF with X number of bins.
## X-HIBF = An HIBF with tmax = X, e.g a maximum of X technical bins on each level.
## ### Column Description ###
## tmax : The maximum number of technical bin on each level
## c_tmax : The technical extra cost of querying an tmax-IBF, compared to 64-IBF
## l_tmax : The estimated query cost for an tmax-HIBF, compared to an 64-HIBF
## m_tmax : The estimated memory consumption for an tmax-HIBF, compared to an 64-HIBF
## (l*m)_tmax : Computed by l_tmax * m_tmax
## size : The expected total size of an tmax-HIBF
# tmax	c_tmax	l_tmax	m_tmax	(l*m)_tmax	size
64	1.00	1.26	1.00	1.26	73.3KiB
128	1.22	1.25	0.66	0.82	48.4KiB
256	1.33	1.33	0.74	0.99	54.3KiB
# Best t_max (regarding expected query runtime): 128
)expected_cout");

    std::string const layout_string{string_from_file(layout_file)};
    EXPECT_NE(layout_string.find("\"tmax\": 128,"), std::string::npos);
}

TEST(execute_estimation_test, with_rearrangement)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const sketches_dir{tmp_dir.path() / "test"};
    std::filesystem::path const input_file{tmp_dir.path() / "test.tsv"};
    std::filesystem::path const layout_file{tmp_dir.path() / "layout.tsv"};
    std::filesystem::path const stats_file{layout_file.string() + ".stats"};

    std::vector<std::string> expected_filenames;
    std::vector<size_t> expected_kmer_counts;

    { // write test.tsv
        std::ofstream fout{input_file};
        for (size_t i{0}; i < 196u;)
        {
            fout << data("seq1.fa").string() << '\t' << (i++) << '\n';
            fout << data("seq2.fa").string() << '\t' << (i++) << '\n';
            fout << data("seq3.fa").string() << '\t' << (i++) << '\n';
            fout << data("small.fa").string() << '\t' << (i++) << '\n';

            expected_filenames.push_back(data("seq1.fa"));
            expected_filenames.push_back(data("seq2.fa"));
            expected_filenames.push_back(data("seq3.fa"));
            expected_filenames.push_back(data("small.fa"));

            expected_kmer_counts.push_back(387);
            expected_kmer_counts.push_back(465);
            expected_kmer_counts.push_back(465);
            expected_kmer_counts.push_back(571);
        }
    }

    chopper::configuration config{};
    config.threads = 1;
    config.k = 15;
    config.column_index_to_cluster = 2;
    config.data_file = input_file;
    config.tmax = 256;
    config.sketch_directory = sketches_dir;
    config.determine_best_tmax = true;
    config.force_all_binnings = true;
    // config.output_verbose_statistics = true;
    config.output_filename = layout_file;

    std::stringstream output_buffer;
    std::stringstream header_buffer;

    chopper::data_store store{.false_positive_rate = config.false_positive_rate,
                              .output_buffer = &output_buffer,
                              .header_buffer = &header_buffer};

    chopper::count::execute(config, store);

    ASSERT_TRUE(std::filesystem::exists(sketches_dir));

    EXPECT_RANGE_EQ(store.filenames, expected_filenames);
    ASSERT_EQ(store.sketches.size(), expected_kmer_counts.size());
    for (size_t i = 0; i < store.sketches.size(); ++i)
    {
        EXPECT_EQ(std::lround(store.sketches[i].estimate()), expected_kmer_counts[i]) << "failed at " << i;

        std::filesystem::path const current_path{expected_filenames[i]};
        std::string const filename = sketches_dir / current_path.stem().string() += ".hll";
        EXPECT_TRUE(std::filesystem::exists(filename));
    }

    chopper::sketch::estimate_kmer_counts(store.sketches, store.kmer_counts);

    chopper::layout::execute(config, store);

    ASSERT_TRUE(std::filesystem::exists(stats_file));

    std::string const written_file{string_from_file(stats_file)};

    EXPECT_EQ(written_file,
              R"expected_cout(## ### Parameters ###
## number of user bins = 196
## number of hash functions = 2
## false positive rate = 0.05
## ### Notation ###
## X-IBF = An IBF with X number of bins.
## X-HIBF = An HIBF with tmax = X, e.g a maximum of X technical bins on each level.
## ### Column Description ###
## tmax : The maximum number of technical bin on each level
## c_tmax : The technical extra cost of querying an tmax-IBF, compared to 64-IBF
## l_tmax : The estimated query cost for an tmax-HIBF, compared to an 64-HIBF
## m_tmax : The estimated memory consumption for an tmax-HIBF, compared to an 64-HIBF
## (l*m)_tmax : Computed by l_tmax * m_tmax
## size : The expected total size of an tmax-HIBF
# tmax	c_tmax	l_tmax	m_tmax	(l*m)_tmax	size
64	1.00	2.23	1.00	2.23	116.5KiB
128	1.22	1.95	1.15	2.24	133.7KiB
256	1.33	1.53	1.19	1.82	138.7KiB
# Best t_max (regarding expected query runtime): 256
)expected_cout");

    std::string const layout_string{string_from_file(layout_file)};
    EXPECT_NE(layout_string.find("\"tmax\": 256,"), std::string::npos);
}
