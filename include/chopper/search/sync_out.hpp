#pragma once

#include <fstream>
#include <mutex>

#include <seqan3/std/filesystem>

class sync_out
{
public:
    sync_out() = default;
    sync_out(sync_out const &) = default;
    sync_out & operator=(sync_out const &) = default;
    sync_out(sync_out &&) = default;
    sync_out & operator=(sync_out &&) = default;
    ~sync_out() = default;

    sync_out(std::filesystem::path const & path) : file(std::ofstream{path}) {}

    template <typename t>
    sync_out & operator<<(t && data)
    {
        std::lock_guard<std::mutex> lock(write_mutex);
        file << std::forward<t>(data);
        return *this;
    }

private:
    std::ofstream file;
    std::mutex write_mutex;
};
