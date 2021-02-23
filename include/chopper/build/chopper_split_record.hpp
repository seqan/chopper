#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include <seqan3/io/sequence_file/all.hpp>

#include <chopper/build/region.hpp>

struct chopper_split_record
{
    std::vector<std::string> filenames{};
    std::vector<size_t> bin_indices{};
    std::vector<size_t> number_of_bins{};
    std::unordered_map<std::string, std::vector<region>> region_map{};

    struct file_type_traits : public seqan3::sequence_file_input_default_traits_dna
    {
        using sequence_alphabet = seqan3::dna4;
    };

    using seq_file_type = seqan3::sequence_file_input<file_type_traits,
                                                      seqan3::fields<seqan3::field::seq, seqan3::field::id>,
                                                      seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

    std::unordered_map<std::string, seqan3::dna4_vector> info;

    void initialise_info() // only initialize info on demand
    {
        for (auto const & filename : filenames)
            for (auto && [seq, id] : seq_file_type{filename})
                info.emplace(filename + id, std::move(seq));
    }

    bool operator==(chopper_split_record const & other) const
    {
        if (std::tie(filenames, bin_indices, number_of_bins) !=
            std::tie(other.filenames, other.bin_indices, other.number_of_bins))
            return false;

        // check region map
        if (region_map.size() != other.region_map.size())
            return false;

        for (auto const & [id, other_regions] : other.region_map)
        {
            auto regions_it = region_map.find(id);
            if (regions_it == region_map.end())
                return false;
            auto const & regions = regions_it->second;

            if (regions.size() != other_regions.size())
                return false;

            for (auto const & reg : other_regions)
                if (std::find(regions.begin(), regions.end(), reg) == regions.end())
                    return false;
        }

        // check info map
        if (info.size() != other.info.size())
            return false;

        for (auto const & [id, other_seq] : other.info)
        {
            auto seq_it = info.find(id);
            if (seq_it == info.end())
                return false;
            auto const & seq = seq_it->second;

            if (seq.size() != other_seq.size())
                return false;

            for (size_t i = 0; i < seq.size(); ++i)
                if (seq[i] != other_seq[i])
                    return false;
        }

        return true;
    }

    bool operator!=(chopper_split_record const & other) const
    {
        return !(*this == other);
    }
};

inline std::ostream & operator<<(std::ostream & s, chopper_split_record const & r)
{
    s << "SPLIT RECORD:"
      << "\n  -> filenames: ";
    for (auto const & filename : r.filenames)
        s << filename << ",";
    s << "\n  -> bin_indices: ";
    for (auto const & index : r.bin_indices)
        s << index << ",";
    s  << "\n  -> number_of_bins: ";
    for (auto const & num : r.number_of_bins)
        s << num << ",";
    s << "\n  -> region_map: ";
    for (auto const & [id, regions] : r.region_map)
    {
        s << '<' << id << ",";
        for (auto const & reg : regions)
            s << reg << '>' << ',';
    }
    s << '\n';

      return s;
}
