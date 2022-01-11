#pragma once

#include <filesystem>

#include <cereal/cereal.hpp>

namespace cereal
{

template <typename archive_t>
void CEREAL_SAVE_FUNCTION_NAME(archive_t & archive, std::filesystem::path const & path)
{
    std::string const str{path.string()};
    archive(str);
}

template <typename archive_t>
void CEREAL_LOAD_FUNCTION_NAME(archive_t & archive, std::filesystem::path & path)
{
    std::string str;
    archive(str);
    path.assign(str);
}

} // namespace cereal
