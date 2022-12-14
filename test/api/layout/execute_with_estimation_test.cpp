#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>

#include "../api_test.hpp"
#include <chopper/count/execute.hpp>
#include <chopper/detail_apply_prefix.hpp>
#include <chopper/layout/execute.hpp>

TEST(execute_estimation_test, few_ubs)
{
    seqan3::test::tmp_filename const input_prefix{"test"};
    seqan3::test::tmp_filename const layout_file{"layout.tsv"};
    std::filesystem::path const stats_file{layout_file.get_path().string() + ".stats"};

    {
        std::ofstream fout{input_prefix.get_path().string() + ".count"};
        fout << "seq0\t500\n"
             << "seq1\t1000\n"
             << "seq2\t500\n"
             << "seq3\t500\n"
             << "seq4\t500\n"
             << "seq5\t500\n"
             << "seq6\t500\n"
             << "seq7\t500\n";
    }

    chopper::configuration config{};
    config.tmax = 64;
    config.determine_best_tmax = true;
    config.input_prefix = input_prefix.get_path();
    config.output_prefix = config.input_prefix;
    config.output_filename = layout_file.get_path();
    chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

    chopper::layout::execute(config);

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
64	1.00	1.00	1.00	1.00	17.5KiB
# Best t_max (regarding expected query runtime): 64
)expected_cout");
}

TEST(execute_estimation_test, many_ubs)
{
    seqan3::test::tmp_filename const input_prefix{"test"};
    seqan3::test::tmp_filename const layout_file{"layout.tsv"};
    std::filesystem::path const stats_file{layout_file.get_path().string() + ".stats"};

    {
        // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
        std::ofstream fout{input_prefix.get_path().string() + ".count"};
        for (size_t i{0}; i < 96u; ++i)
            fout << seqan3::detail::to_string("seq", i, '\t', 100 * ((i + 20) / 20), '\n');
    }

    chopper::configuration config{};
    config.tmax = 1024;
    config.determine_best_tmax = true;
    config.input_prefix = input_prefix.get_path();
    config.output_prefix = config.input_prefix;
    config.output_filename = layout_file.get_path();
    chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

    chopper::layout::execute(config);

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
64	1.00	1.26	1.00	1.26	85.2KiB
128	0.96	0.98	0.59	0.58	50.3KiB
256	1.20	1.20	0.70	0.83	59.3KiB
# Best t_max (regarding expected query runtime): 128
)expected_cout");

    std::string const layout_string{string_from_file(layout_file.get_path())};
    EXPECT_NE(layout_string.find("\"tmax\": 128,"), std::string::npos);
}

TEST(execute_estimation_test, many_ubs_force_all)
{
    seqan3::test::tmp_filename const input_prefix{"test"};
    seqan3::test::tmp_filename const layout_file{"layout.tsv"};
    std::filesystem::path const stats_file{layout_file.get_path().string() + ".stats"};

    {
        // There are 20 files with a count of {100,200,300,400} each. There are 16 files with count 500.
        std::ofstream fout{input_prefix.get_path().string() + ".count"};
        for (size_t i{0}; i < 96u; ++i)
            fout << seqan3::detail::to_string("seq", i, '\t', 100 * ((i + 20) / 20), '\n');
    }

    chopper::configuration config{};
    config.tmax = 256;
    config.determine_best_tmax = true;
    config.force_all_binnings = true;
    config.input_prefix = input_prefix.get_path();
    config.output_prefix = config.input_prefix;
    config.output_filename = layout_file.get_path();
    chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

    chopper::layout::execute(config);

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
64	1.00	1.26	1.00	1.26	85.2KiB
128	0.96	0.98	0.59	0.58	50.3KiB
256	1.20	1.20	0.70	0.83	59.3KiB
# Best t_max (regarding expected query runtime): 128
)expected_cout");

    std::string const layout_string{string_from_file(layout_file.get_path())};
    EXPECT_NE(layout_string.find("\"tmax\": 128,"), std::string::npos);
}

TEST(execute_estimation_test, with_rearrangement)
{
    seqan3::test::tmp_filename const prefix{"test"};
    seqan3::test::tmp_filename const input_file{"test.tsv"};
    seqan3::test::tmp_filename const layout_file{"layout.tsv"};
    std::filesystem::path const stats_file{layout_file.get_path().string() + ".stats"};

    {
        std::ofstream fout{input_file.get_path()};
        for (size_t i{0}; i < 196u;)
        {
            fout << data("seq1.fa").string() << '\t' << (i++) << '\n';
            fout << data("seq2.fa").string() << '\t' << (i++) << '\n';
            fout << data("seq3.fa").string() << '\t' << (i++) << '\n';
            fout << data("small.fa").string() << '\t' << (i++) << '\n';
        }
    }

    {
        chopper::configuration config{};
        config.threads = 1;
        config.k = 15;
        config.column_index_to_cluster = 2;
        config.data_file = input_file.get_path();
        config.output_prefix = prefix.get_path();
        chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

        chopper::count::execute(config);
    }

    ASSERT_TRUE(std::filesystem::exists(prefix.get_path().string() + ".count"));
    ASSERT_TRUE(std::filesystem::exists(prefix.get_path().string() + "_sketches"));
    ASSERT_TRUE(std::filesystem::exists(prefix.get_path().string() + "_sketches/seq1.hll"));
    ASSERT_TRUE(std::filesystem::exists(prefix.get_path().string() + "_sketches/seq2.hll"));
    ASSERT_TRUE(std::filesystem::exists(prefix.get_path().string() + "_sketches/seq3.hll"));
    ASSERT_TRUE(std::filesystem::exists(prefix.get_path().string() + "_sketches/small.hll"));

    chopper::configuration config{};
    config.tmax = 256;
    config.rearrange_user_bins = true;
    config.determine_best_tmax = true;
    config.force_all_binnings = true;
    // config.output_verbose_statistics = true;
    config.input_prefix = prefix.get_path();
    config.output_prefix = config.input_prefix;
    config.output_filename = layout_file.get_path();
    chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

    chopper::layout::execute(config);

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
128	0.96	1.57	1.15	1.81	134.3KiB
256	1.20	1.39	1.19	1.66	138.7KiB
# Best t_max (regarding expected query runtime): 256
)expected_cout");

    std::string const layout_string{string_from_file(layout_file.get_path())};
    EXPECT_NE(layout_string.find("\"tmax\": 256,"), std::string::npos);
}
