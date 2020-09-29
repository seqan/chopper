#pragma once

#include <seqan3/core/debug_stream.hpp>

// helper function to print a matrix when debugging
template <typename matrix_type,typename matrix_value_type>
void print_matrix(matrix_type const & matrix,
                  size_t const row_bound,
                  size_t const column_bound,
                  matrix_value_type const inf)
{
    for (size_t i = 0; i < row_bound; ++i)
    {
        for (size_t j = 0; j < column_bound; ++j)
        {
            if (matrix[i][j] == inf)
                seqan3::debug_stream << "inf\t";
            else
                seqan3::debug_stream << matrix[i][j] << '\t';
        }
        seqan3::debug_stream << '\n';
    }
    seqan3::debug_stream << '\n';
}
