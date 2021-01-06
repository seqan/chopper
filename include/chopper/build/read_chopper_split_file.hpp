#pragma once

#include <cstdlib>
#include <fstream>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <seqan3/std/charconv>

#include <chopper/build/batch.hpp>
#include <chopper/build/build_data.hpp>
#include <chopper/build/region.hpp>
#include <chopper/detail_bin_prefixes.hpp>
#include <chopper/detail_starts_with.hpp>

void parse_chopper_split_header_line(std::string const & line, build_data & data)
{
    if (line.substr(1, hibf_prefix.size()) == hibf_prefix)
    {
        assert(line.substr(hibf_prefix.size() + 2, 11) == "max_bin_id:");
        auto it = std::find(line.begin() + hibf_prefix.size() + 13, line.end(), '_'); // skip "MERGED"/"SPLIT"
        ++it; // skip "_"
        it = std::find(it, line.end(), '_'); // skip "BIN"
        ++it; // skip "_"
        data.hibf_max_bin = std::atoi(std::string(it, line.end()).c_str());
    }
    else if (line.substr(1, merged_bin_prefix.size()) == merged_bin_prefix)
    {
        std::string const hidx_str(line.begin() + 1 /*#*/ + merged_bin_prefix.size() + 1 /*_*/,
                                   std::find(line.begin() + merged_bin_prefix.size() + 2, line.end(), ' '));
        assert(line.substr(merged_bin_prefix.size() + hidx_str.size() + 3, 11) == "max_bin_id:");
        std::string const lidx_str = line.substr(merged_bin_prefix.size() + hidx_str.size() + 14,
                                                 line.size() - merged_bin_prefix.size() - hidx_str.size() - 14);

        size_t const hidx = std::atoi(hidx_str.c_str());
        size_t const lidx = std::atoi(lidx_str.c_str());

        data.merged_max_bin_map.emplace(hidx, lidx);
    }
}

auto parse_chopper_split_line(std::string const & line)
{
    char const * buffer = line.c_str();
    auto field_start = &buffer[0];
    auto field_end = &buffer[0];
    auto const buffer_end = field_start + line.size();

    while (field_end != buffer_end && *field_end != '\t') ++field_end;

    std::string const filename(field_start, field_end);

    assert(*field_end == '\t');
    ++field_end; // skip tab
    field_start = field_end;
    while (field_end != buffer_end && *field_end != '\t') ++field_end;

    std::string const id(field_start, field_end);

    region reg{};

    // read begin
    assert(*field_end == '\t');
    ++field_end;
    auto res = std::from_chars(field_end, buffer_end, reg.begin);
    field_end = res.ptr;

    // read end
    assert(*field_end == '\t');
    ++field_end;
    res = std::from_chars(field_end, buffer_end, reg.end);
    field_end = res.ptr;

    // read hibf bin index
    assert(*field_end == '\t');
    ++field_end;
    res = std::from_chars(field_end, buffer_end, reg.hidx);
    field_end = res.ptr;

    // read libf bin index
    assert(*field_end == '\t');
    ++field_end;

    if (*field_end == '-')
        reg.lidx = -1;
    else
        res = std::from_chars(field_end, buffer_end, reg.lidx);

    return std::make_tuple(std::move(filename), std::move(id), std::move(reg));
}

auto read_chopper_split_file(std::string const & chopper_split_filename)
{
    build_data data;
    std::vector<batch> batches{};

    std::ifstream chopper_split_file{chopper_split_filename};

    if (!chopper_split_file.good() || !chopper_split_file.is_open())
        throw std::logic_error{"Could not open file for reading"};

    // parse data
    // -------------------------------------------------------------------------
    std::string current_line;
    while (std::getline(chopper_split_file, current_line) && current_line[0] == '#')
        parse_chopper_split_header_line(current_line, data);

    // parse lines into vectors which are parsed into batches in the next step
    // -------------------------------------------------------------------------
    std::vector<std::set<std::string>> files_per_hbin{}; // temporary storage for filenames per high level ibf bin
    std::vector<int64_t> libf_bins_per_hbin{}; // temporary storage for number of low level bins if any
    do
    {
        auto const && [filename, seq_id, region_record] = parse_chopper_split_line(current_line);
        data.num_libfs += (region_record.lidx != -1);

        std::string const file_seq_id{filename + seq_id};
        data.region_map[file_seq_id].push_back(region_record);

        if (files_per_hbin.size() <= region_record.hidx)
            files_per_hbin.resize(region_record.hidx + 1);

        if (libf_bins_per_hbin.size() <= region_record.hidx)
            libf_bins_per_hbin.resize(region_record.hidx + 1);

        files_per_hbin[region_record.hidx].insert(filename);
        libf_bins_per_hbin[region_record.hidx] = std::max(libf_bins_per_hbin[region_record.hidx], region_record.lidx + 1);

    } while (std::getline(chopper_split_file, current_line));

    // process files_per_hbin into batches
    // -------------------------------------------------------------------------
    // Note: The following code is probably not the most efficient one, but it is reliable
    assert(files_per_hbin.size() != 0); // there should be some bins in the HIBF
    int64_t hibf_max_batch_record_pos{-1};

    // 1. Each merged bin is a batch
    for (size_t hidx = 0; hidx < files_per_hbin.size(); ++hidx)
    {
        if (libf_bins_per_hbin[hidx] > 0) // there are LIBF bins. This indicates a merged bin
        {
            if (hidx == data.hibf_max_bin)
                hibf_max_batch_record_pos = batches.size();

            std::vector<std::string> filenames(files_per_hbin[hidx].begin(), files_per_hbin[hidx].end());
            assert(libf_bins_per_hbin[hidx] >= 0);
            size_t num_libfs = static_cast<size_t>(libf_bins_per_hbin[hidx]);
            batches.push_back(batch{std::move(filenames), {hidx}, num_libfs});

            files_per_hbin[hidx].clear(); // empty out bin, s.t. it doesn't participate below
        }
    }

    // 2. Files can be associated with multiple bins and bins can be associated with multiple files.
    //    Therefore all files and bins that are connected somehow would go into one batch to avoid,
    //    that a File and it's sequences are processed more than once.
    for (size_t hidx = 0; hidx < files_per_hbin.size(); ++hidx)
    {
        if (files_per_hbin[hidx].empty())
            continue;

        std::set<std::string> & current_filename_set{files_per_hbin[hidx]};
        std::vector<size_t> current_hidxs{hidx};

        assert(libf_bins_per_hbin[hidx] == 0); // merged bins should have been tackled before

        for (size_t i = hidx + 1; i < files_per_hbin.size(); ++i)
        {
            for (auto & filename : files_per_hbin[i])
            {
                if (current_filename_set.find(filename) != current_filename_set.end())
                {
                    current_filename_set.merge(files_per_hbin[i]);
                    current_hidxs.push_back(i);
                    files_per_hbin[i].clear();
                    break;
                }
            }
        }

        if (std::find(current_hidxs.begin(), current_hidxs.end(), data.hibf_max_bin) != current_hidxs.end())
            hibf_max_batch_record_pos = batches.size();

        std::vector<std::string> filenames(current_filename_set.begin(), current_filename_set.end());
        batches.push_back(batch{std::move(filenames), std::move(current_hidxs), 0u});

    }

    assert(hibf_max_batch_record_pos != -1);
    data.hibf_max_batch_record = &batches[hibf_max_batch_record_pos];
    data.hibf_num_technical_bins = files_per_hbin.size();
    return std::make_pair(std::move(data), std::move(batches));
};
