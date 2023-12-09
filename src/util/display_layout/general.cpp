// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/input.hpp>

#include <hibf/contrib/robin_hood.hpp>
#include <hibf/contrib/std/chunk_by_view.hpp>
#include <hibf/contrib/std/to.hpp>
#include <hibf/sketch/hyperloglog.hpp>

#include "shared.hpp"

struct progress_bar
{
    progress_bar(size_t const total) : total{total}
    {
        print_progress(0);
    }

    void report()
    {
        // Locks on construction, unlocks on deconstruction. RAII-way to handle mutex lock/unlock.
        // Only the thread which acquired the lock can proceed. The other threads have to **wait** until
        // the lock is released.
        std::lock_guard report_guard{report_mutex};
        ++current;

        // Would also work fine without this special case. However, we only want to announce "100% finished" once we
        // are actually finished.
        if (current == total)
        {
            print_progress(100);
            return;
        }

        size_t const percentage = std::min<size_t>(99, std::ceil(100 * current / static_cast<double>(total)));
        if (percentage > last_printed)
        {
            last_printed = percentage;
            print_progress(percentage);
        }
    }

    void print_progress(size_t const percentage) const
    {
        assert(percentage <= 100u);
        std::cerr << '[';
        for (size_t i{}; i < percentage; ++i)
            std::cerr << '=';
        if (percentage < 100u)
        {
            std::cerr << '>';
            for (size_t i{1u}; i < 100u - percentage; ++i)
                std::cerr << ' ';
        }
        std::cerr << "] " << percentage << " %\r" << std::flush;
    }

    std::mutex report_mutex{};
    size_t last_printed{};
    size_t current{};
    size_t total{};
};

// Using two sets and erasing from shared is slower
void keep_duplicates(robin_hood::unordered_set<uint64_t> & shared, std::vector<uint64_t> const & current)
{
    // static + calling clear is slower
    robin_hood::unordered_set<uint64_t> result{};

    for (uint64_t value : current)
    {
        if (shared.contains(value))
            result.emplace(value);
    }

    shared = std::move(result);
}

struct record
{
    size_t tb_index{};
    size_t exact_size{};
    size_t estimated_size{};
    size_t shared_size{};
    size_t ub_count{};
    std::string_view kind{};
    size_t splits{};

    void write_to(std::ostream & stream) const
    {
        size_t split_count{splits};
        for (size_t i{}, total{split_count}; i < total; ++i)
        {
            stream << tb_index + i << '\t'   //
                   << exact_size << '\t'     //
                   << estimated_size << '\t' //
                   << shared_size << '\t'    //
                   << ub_count << '\t'       //
                   << kind << '\t'           //
                   << split_count << '\n';
            // Subsequent split bins display 0, the first split bin displays the actual split count.
            split_count = 0u;
        }
    }

    static void write_header_to(std::ostream & stream, std::string_view const layout_filename)
    {
        stream << "# Layout: " << layout_filename << '\n' //
               << "tb_index\t"
               << "exact_size\t"
               << "estimated_size\t"
               << "shared_size\t"
               << "ub_count\t"
               << "kind\t"
               << "splits\n";
        stream << std::flush;
    }
};

void process_and_write_records_to(std::vector<record> & records, std::ostream & stream)
{
    std::ranges::sort(records,
                      [](record const & lhs, record const & rhs)
                      {
                          return lhs.tb_index < rhs.tb_index;
                      });
    // Split bins offset the tb index. At this point of code, the tb index is just consecutive and does not
    // count split bins:
    // tb_index ub_count kind   split_count
    // 0        4        merged 1
    // 1        1        split  3
    // 2        3        merged 1
    // Without the correction, the output would be:
    // 0        4        merged 1
    // 1        1        split  3
    // 2        1        split  0
    // 3        1        split  0
    // 2        3        merged 1
    // Actually, the split bin will occupy 3 technical bins, i.e. we want:
    // 0        4        merged 1
    // 1        1        split  3
    // 2        1        split  0
    // 3        1        split  0
    // 4        3        merged 1
    size_t tb_offset{};
    for (auto & record : records)
    {
        assert(record.splits > 0u);
        record.tb_index += tb_offset;
        tb_offset += record.splits - 1u;
    }
    // Now we can print the results.
    for (auto & record : records)
        record.write_to(stream);

    stream << std::flush;
}

int execute(config const & cfg)
{
    std::ifstream layout_file{cfg.input};

    if (!layout_file.good() || !layout_file.is_open())
        throw std::logic_error{"Could not open file " + cfg.input.string() + " for reading"};

// https://godbolt.org/z/PeKnxzjn1
#if defined(__clang__)
    auto tuple = chopper::layout::read_layouts_file(layout_file);
    // https://godbolt.org/z/WoWf55KPb
    auto filenames = std::move(std::get<0>(tuple));
    auto chopper_config = std::move(std::get<1>(tuple));
    auto hibf_layouts = std::move(std::get<2>(tuple));
#else
    auto [filenames, chopper_config, hibf_layouts] = chopper::layout::read_layouts_file(layout_file);
#endif
    auto const & hibf_config = chopper_config.hibf_config;

    layout_file.close();
    std::ofstream output_stream{cfg.output};

    if (!output_stream.good() || !output_stream.is_open())
        throw std::logic_error{"Could not open file " + cfg.output.string() + " for reading"};

    record::write_header_to(output_stream, cfg.input.c_str());

    // Fetch all file sizes such that sorting by file size doesn't have to access the filesystem too often.
    // n = filenames.size()
    // Constructing this vector has `n` filesystem accesses.
    // Sorting without pre-fetching has `O(n * log(n))` accesses.
    std::vector<std::uintmax_t> const filesizes{[&filenames]()
                                                {
                                                    std::vector<std::uintmax_t> result{};
                                                    result.reserve(filenames.size());
                                                    for (auto const & filename : filenames)
                                                        result.push_back(std::filesystem::file_size(filename.front()));
                                                    return result;
                                                }()};

    // Sorts by the technical bin indices in the top-level IBF:
    // split bins: storage_TB_id, previous_TB_indices is empty
    // merged bins: previous_TB_indices[0]
    // If the index is the same, sort by file sizes (happens for merged bins).
    // Using the smallest file to initialise the shared k-mers later will be less work.
    std::ranges::sort(
        hibf_layouts[0].user_bins,
        [&filesizes](seqan::hibf::layout::layout::user_bin const & lhs,
                     seqan::hibf::layout::layout::user_bin const & rhs)
        {
            size_t const first_idx = lhs.previous_TB_indices.empty() ? lhs.storage_TB_id : lhs.previous_TB_indices[0];
            size_t const second_idx = rhs.previous_TB_indices.empty() ? rhs.storage_TB_id : rhs.previous_TB_indices[0];
            return first_idx < second_idx || (first_idx == second_idx && filesizes[lhs.idx] < filesizes[rhs.idx]);
        });

    size_t const total_ub_count = hibf_layouts[0].user_bins.size();
    progress_bar progress{total_ub_count};

    // Create chunks containing user bin indices for one technical bin.
    // We need to convert the chunk_by_view to a std::vector because we need a random access range for the
    // parallelisation. As a side effect, std::vector is also a sized range, which makes things even easier.
    // E.g., `[0,1,2,3] [4] [5,6,7,8]`:
    // TB UBs
    // 0  0,1,2,3
    // 1  4
    // 2  5,6,7,8
    std::vector<std::vector<size_t>> const chunks = [&]()
    {
        auto ub_indices = std::views::iota(size_t{}, total_ub_count);
        // Two user bins belong to the same chunk if they are in the same technical bin.
        auto predicate = [&](size_t const lhs, size_t const rhs)
        {
            auto const & lhs_ub = hibf_layouts[0].user_bins[lhs];
            auto const & rhs_ub = hibf_layouts[0].user_bins[rhs];
            // The top-level technical bin index for the current user bin.
            // user_bin.previous_TB_indices.size() == 0: true for split bins, false for merged bins
            // user_bin.storage_TB_id: technical bin index on the lowest level
            // user_bin.previous_TB_indices[0]: technical bin index on the top-level
            size_t const lhs_idx =
                (lhs_ub.previous_TB_indices.size() == 0) ? lhs_ub.storage_TB_id : lhs_ub.previous_TB_indices[0];
            size_t const rhs_idx =
                (rhs_ub.previous_TB_indices.size() == 0) ? rhs_ub.storage_TB_id : rhs_ub.previous_TB_indices[0];

            return lhs_idx == rhs_idx;
        };
        auto chunked_by_tb = ub_indices | seqan::stl::views::chunk_by(std::move(predicate));

        return seqan::stl::ranges::to<std::vector<std::vector<size_t>>>(chunked_by_tb);
    }();

    std::vector<record> records(chunks.size());

#pragma omp parallel for schedule(dynamic) num_threads(cfg.threads)
    for (size_t tb_index = 0; tb_index < chunks.size(); ++tb_index)
    {
        auto const & chunk = chunks[tb_index];
        // Estimates the cardinality of one technical bin. For merged bins, user bins will be iteratively added.
        seqan::hibf::sketch::hyperloglog sketch{hibf_config.sketch_bits};
        // Used to determine the exact cardinality for one technical bin.
        robin_hood::unordered_set<uint64_t> current_kmer_set{};
        // Stores shared k-mers across user bins of a merged technical bin.
        robin_hood::unordered_set<uint64_t> shared_kmers{};
        // We can't use `shared_kmers.size() == 0` instead of `shared_kmers_initialised`, because keep_duplicates
        // will result in a size of 0 when there are no shared k-mers.
        bool shared_kmers_initialised{false};
        std::vector<uint64_t> current_kmers{};
        // How many user bins are stored in the current technical bin? Always 1 for split bins.
        size_t const ub_count{chunk.size()};
        bool const is_merged{ub_count > 1u};

        for (size_t const ub_index : chunk)
        {
            auto const & user_bin = hibf_layouts[0].user_bins[ub_index];
            current_kmers.clear();

            // We don't need to keep the current_kmers if there are no shared k-mers to merge them with.
            bool const fill_current_kmers = is_merged && !(shared_kmers_initialised && shared_kmers.empty());

            for (auto const & filename : filenames[user_bin.idx])
            {
                process_file(filename,
                             current_kmer_set,
                             current_kmers,
                             sketch,
                             fill_current_kmers,
                             chopper_config.k,
                             chopper_config.window_size);
            }

            // Compute set intersection: shared_kmers = shared_kmers ∩ current_kmers
            // This happens for each user bin that belongs to a merged bin.
            if (fill_current_kmers)
            {
                if (!shared_kmers_initialised)
                {
                    shared_kmers_initialised = true;
                    shared_kmers.insert(current_kmers.begin(), current_kmers.end());
                }
                else
                {
                    keep_duplicates(shared_kmers, current_kmers);
                }
            }

            progress.report();
        }

        // Into how many techincal bins is the user bin split? Always 1 for merged bins.
        size_t const split_count{is_merged ? 1u : hibf_layouts[0].user_bins[chunk[0]].number_of_technical_bins};
        size_t const avg_kmer_count = (current_kmer_set.size() + split_count - 1u) / split_count;
        size_t const sketch_estimate = (sketch.estimate() + split_count - 1u) / split_count;

        records[tb_index] = record{.tb_index = tb_index,
                                   .exact_size = avg_kmer_count,
                                   .estimated_size = sketch_estimate,
                                   .shared_size = shared_kmers.size(),
                                   .ub_count = ub_count,
                                   .kind = (is_merged ? "merged" : "split"),
                                   .splits = split_count};
    }

    process_and_write_records_to(records, output_stream);

    // The progress bar uses a carriage return '\r' to only use a single line.
    std::cerr << '\n';

    return 0;
}

void execute_general(config const & cfg)
{
    execute(cfg);
}
