#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <chopper/layout/execute.hpp>
#include <chopper/sketch/estimate_kmer_counts.hpp>
#include <chopper/sketch/execute.hpp>

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

    chopper::layout::layout hibf_layout{};

    chopper::data_store store{.false_positive_rate = config.false_positive_rate,
                              .hibf_layout = &hibf_layout,
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

    chopper::layout::layout hibf_layout{};

    chopper::data_store data{.false_positive_rate = config.false_positive_rate,
                             .hibf_layout = &hibf_layout,
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
    EXPECT_EQ(layout_string,
              R"expected_layout(##CONFIG:
##{
##    "config": {
##        "version": 2,
##        "data_file": {
##            "value0": ""
##        },
##        "debug": false,
##        "sketch_directory": {
##            "value0": ""
##        },
##        "k": 19,
##        "sketch_bits": 12,
##        "disable_sketch_output": true,
##        "precomputed_files": false,
##        "output_filename": {
##            "value0": ")expected_layout"
                  + layout_file.string() +
                  R"expected_layout("
##        },
##        "tmax": 128,
##        "num_hash_functions": 2,
##        "false_positive_rate": 0.05,
##        "alpha": 1.2,
##        "max_rearrangement_ratio": 0.5,
##        "threads": 1,
##        "disable_estimate_union": true,
##        "disable_rearrangement": true,
##        "determine_best_tmax": true,
##        "force_all_binnings": false
##    }
##}
##ENDCONFIG
#HIGH_LEVEL_IBF max_bin_id:14
#MERGED_BIN_14 max_bin_id:0
#MERGED_BIN_15 max_bin_id:0
#FILES	BIN_INDICES	NUMBER_OF_BINS
seq0	0	1
seq19	1	1
seq18	2	1
seq17	3	1
seq16	4	1
seq15	5	1
seq14	6	1
seq13	7	1
seq12	8	1
seq11	9	1
seq10	10	1
seq9	11	1
seq8	12	1
seq7	13	1
seq4	14;0	1;22
seq5	14;22	1;21
seq6	14;43	1;21
seq1	15;0	1;22
seq2	15;22	1;21
seq3	15;43	1;21
seq32	16	1
seq33	17	1
seq34	18	1
seq35	19	1
seq30	20	1
seq36	21	1
seq37	22	1
seq38	23	1
seq39	24	1
seq31	25	1
seq20	26	1
seq21	27	1
seq22	28	1
seq23	29	1
seq25	30	1
seq26	31	1
seq27	32	1
seq28	33	1
seq24	34	1
seq29	35	1
seq40	36	1
seq41	37	1
seq42	38	1
seq43	39	1
seq44	40	1
seq45	41	1
seq46	42	1
seq47	43	1
seq48	44	1
seq59	45	1
seq58	46	1
seq57	47	1
seq56	48	1
seq55	49	1
seq54	50	1
seq53	51	1
seq52	52	1
seq51	53	1
seq50	54	1
seq49	55	1
seq79	56	2
seq78	58	2
seq77	60	2
seq76	62	2
seq75	64	2
seq74	66	2
seq73	68	2
seq70	70	2
seq71	72	2
seq72	74	2
seq60	76	2
seq61	78	2
seq62	80	2
seq63	82	2
seq64	84	2
seq65	86	2
seq66	88	2
seq67	90	2
seq68	92	2
seq69	94	2
seq88	96	2
seq94	98	2
seq93	100	2
seq92	102	2
seq91	104	2
seq90	106	2
seq89	108	2
seq87	110	2
seq86	112	2
seq85	114	2
seq84	116	2
seq83	118	2
seq82	120	2
seq81	122	2
seq80	124	2
seq95	126	2
)expected_layout");
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

    chopper::layout::layout hibf_layout{};

    chopper::data_store data{.false_positive_rate = config.false_positive_rate,
                             .hibf_layout = &hibf_layout,
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

    chopper::layout::layout hibf_layout{};

    chopper::data_store store{.false_positive_rate = config.false_positive_rate, .hibf_layout = &hibf_layout};

    chopper::sketch::execute(config, store);

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
