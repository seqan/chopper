#pragma once

#include <string>
#include <vector>

struct data_file_record
{
    /*!\name Construction, destruction and assignment
     * \{
     */
    data_file_record() = default;
    data_file_record(data_file_record const & rhs) = default;
    data_file_record(data_file_record && rhs) = default;
    data_file_record & operator=(data_file_record const & rhs) = default;
    data_file_record & operator=(data_file_record && rhs) = default;
    ~data_file_record() = default;

    data_file_record(std::string n, std::vector<std::string> f, size_t b) :
        bin_name{std::move(n)},
        filenames{std::move(f)},
        bins{b}
    {}
    //!\}

    std::string bin_name{};
    std::vector<std::string> filenames{};
    size_t bins{};
};
