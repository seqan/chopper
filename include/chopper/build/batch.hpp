#pragma once

#include <string>
#include <vector>

struct batch
{
    std::vector<std::string> filenames{};
    std::vector<size_t> hibf_bins{};
    size_t libf_num_bins{};

    bool operator==(batch const & b) const
    {
        return std::tie(filenames, hibf_bins, libf_num_bins) == std::tie(b.filenames, b.hibf_bins, b.libf_num_bins);
    }

    bool operator!=(batch const & b) const
    {
        return std::tie(filenames, hibf_bins, libf_num_bins) != std::tie(b.filenames, b.hibf_bins, b.libf_num_bins);
    }
};

inline std::ostream & operator<<(std::ostream & s, batch const & b)
{
    s << "BATCH: \n"
      << "  -> filenames: ";
    for (auto const & filename : b.filenames)
        s << filename << ",";
    s << "\n  -> hidxs: ";
    for (auto const & hidx : b.hibf_bins)
        s << hidx << ",";
    s << "\n  -> libf_num_bins: " << b.libf_num_bins << '\n';

    return s;
}
