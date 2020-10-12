#pragma once

#include <seqan/graph_msa.h>

// Helper function for rounding to n significant digits.  Ported from
// Java code found here: http://stackoverflow.com/questions/202302
inline double roundToSignificantFigures(double num, int n)
{
    if (num == 0)
        return 0;

    const double d = ceil(log10(num < 0 ? -num : num));
    const int power = n - (int) d;

    const double magnitude = pow(10.0, power);
    const long shifted = static_cast<long>(round(num*magnitude));
    return shifted / magnitude;
}

template<typename TMatrix>
auto neighbour_joining(TMatrix mat)
{
    using TSize = typename seqan::Size<TMatrix>::Type;
    using node_weight_type = double;
    using graph_type = seqan::Graph<seqan::Tree<node_weight_type>>;
    using TVertexDescriptor = typename seqan::VertexDescriptor<graph_type>::Type;

    graph_type g; // resulting guide tree

    TVertexDescriptor nilVertex = seqan::getNil<TVertexDescriptor>();
    TSize nseq = (TSize) std::sqrt((double)length(mat));

    // Assert that the input matrix has no negative values.
// #if SEQAN_ENABLE_DEBUG
//     for (unsigned i = 0; i < length(mat); ++i)
//         SEQAN_ASSERT_GEQ_MSG(mat[i], 0, "i = %u", i);
// #endif  // #if SEQAN_ENABLE_DEBUG

    //for(TSize i=0;i<nseq;++i) {
    //    for(TSize j=0;j<nseq;++j) {
    //        std::cout << seqan::getValue(mat, i*nseq+j) << ",";
    //    }
    //    std::cout << std::endl;
    //}

    // Handle base cases for one and two sequences.
    seqan::clearVertices(g);
    if (nseq == 1)
    {
        g.data_root = seqan::addVertex(g);
        return g;
    }
    else if (nseq == 2)
    {
        TVertexDescriptor v1 = seqan::addVertex(g);
        TVertexDescriptor v2 = seqan::addVertex(g);
        TVertexDescriptor internalVertex = seqan::addVertex(g);
        seqan::addEdge(g, internalVertex, v1, (node_weight_type) roundToSignificantFigures(mat[1] / 2.0, 5));
        seqan::addEdge(g, internalVertex, v2, (node_weight_type) roundToSignificantFigures(mat[1] / 2.0, 5));
        g.data_root = internalVertex;
        return g;
    }

    // Create a normalized copy of mat with fixed point numbers, precision of 10 digits.
    double normFactor = 1.0;
    for (unsigned i = 0; i < length(mat); ++i)
        mat[i] = static_cast<int64_t>(10000000.0 * ((double)(mat[i]) / (double)(normFactor)));

    // First initialization
    seqan::String<int64_t> av;    // Average branch length to a combined node
    seqan::resize(av, nseq, 0);

    seqan::String<TVertexDescriptor> connector;   // Nodes that need to be connected
    seqan::resize(connector, nseq);

    for (TSize i = 0; i < nseq; ++i)
    {
        seqan::addVertex(g);  // Add all the nodes that correspond to sequences
        connector[i] = i;
        mat[i*nseq+i] = 0;
    }

    int64_t sumOfBranches = 0;
    // Determine the sum of all branches and
    // copy upper triangle mat to lower triangle
    for (TSize col = 1; col < nseq; ++col)
        for (TSize row = 0; row < col; ++row)
            sumOfBranches += mat[col*nseq+row] = mat[row*nseq+col];

    // Main cycle
    int64_t fnseqs = static_cast<int64_t>(nseq);
    for (TSize nc = 0; nc < (nseq - 3); ++nc)
    {
        // Compute the sum of branch lengths for all possible pairs
        bool notFound = true;
        int64_t tmin = 0;
        TSize mini = 0;  // Next pair of seq i and j to join
        TSize minj = 0;
        int64_t diToAllOthers = 0;
        int64_t djToAllOthers = 0;
        int64_t total = 0;
        int64_t dMinIToOthers = 0;
        int64_t dMinJToOthers = 0;
        for (TSize col = 1; col < nseq; ++col)
        {
            if (connector[col] != nilVertex)
            {
                for (TSize row = 0; row < col; ++row)
                {
                    if (connector[row] != nilVertex)
                    {
                        diToAllOthers = 0;
                        djToAllOthers = 0;

                        for (TSize i = 0; i < nseq; ++i)
                        {
                            diToAllOthers += mat[i*nseq+row];
                            djToAllOthers += mat[i*nseq+col];
                        }

                        total = diToAllOthers + djToAllOthers + (fnseqs - 2) * mat[row*nseq+col] + 2 * (sumOfBranches - diToAllOthers - djToAllOthers);
                        total /= (2*(fnseqs - 2));

                        if ((notFound) || (total < tmin))
                        {
                            notFound = false;
                            tmin = total;
                            mini = row;
                            minj = col;
                            dMinIToOthers = diToAllOthers;
                            dMinJToOthers = djToAllOthers;
                        }
                    }
                }
            }
        }

        // Print nodes that are about to be joined
        //std::cout << mini << std::endl;
        //std::cout << minj << std::endl;
        //std::cout << tmin << std::endl;
        //std::cout << std::endl;

        int64_t dmin = mat[mini*nseq + minj];
        dMinIToOthers = dMinIToOthers / (fnseqs - 2);
        dMinJToOthers = dMinJToOthers / (fnseqs - 2);
        int64_t iBranch = (dmin + dMinIToOthers - dMinJToOthers) / 2;
        int64_t jBranch = dmin - iBranch;
        iBranch -= av[mini];
        jBranch -= av[minj];

        // Set negative branch length to zero
        if (iBranch < 0) iBranch = 0;
        if (jBranch < 0) jBranch = 0;

        // Print branch lengths
        //std::cout << iBranch << std::endl;
        //std::cout << jBranch << std::endl;
        //std::cout << std::endl;

        // Build tree
        TVertexDescriptor internalVertex = seqan::addVertex(g);
        seqan::addEdge(g, internalVertex, connector[mini], (node_weight_type) roundToSignificantFigures((iBranch / 10000000.0) * normFactor, 5));
        seqan::addEdge(g, internalVertex, connector[minj], (node_weight_type) roundToSignificantFigures((jBranch / 10000000.0) * normFactor, 5));

        // Remember the average branch length for the new combined node
        // Must be subtracted from all branches that include this node
        if (dmin < 0)
            dmin = 0;
        av[mini] = dmin / 2;

        // Re-initialisation
        // mini becomes the new combined node, minj is killed
        --fnseqs;
        connector[minj] = nilVertex;
        connector[mini] = internalVertex;

        for (TSize j = 0; j < nseq; ++j)
        {
            if (connector[j] != nilVertex && mini != j)
            {
                TSize const new_value = (mat[mini*nseq+j] + mat[minj*nseq+j]) / 2;
                assert(mat[mini*nseq+j] == mat[j*nseq+mini]);
                mat[mini*nseq+j] = new_value;
                mat[j*nseq+mini] = new_value;
                sumOfBranches -= new_value;
            }
            else
            {
                assert(mat[j*nseq+minj] == mat[minj*nseq+j]);
                sumOfBranches -= mat[j*nseq+minj];
            }

            mat[j*nseq+minj] = mat[minj*nseq+j] = 0;
        }
    }

    // Only three nodes left

    // Find the remaining nodes
    seqan::String<TSize> l;
    seqan::resize(l,3);
    TSize count = 0;
    for (TSize i=0; i<nseq; ++i)
    {
        if (connector[i] != nilVertex)
        {
            l[count] = i;
            ++count;
        }
    }

    // Remaining nodes
    //std::cout << l[0] << std::endl;
    //std::cout << l[1] << std::endl;
    //std::cout << l[2] << std::endl;
    //std::cout << std::endl;

    seqan::String<int64_t> branch;
    seqan::resize(branch, 3);
    branch[0] = (mat[l[0]*nseq+l[1]] + mat[l[0]*nseq+l[2]] - mat[l[1]*nseq+l[2]]) / 2;
    branch[1] = (mat[l[1]*nseq+l[2]] + mat[l[0]*nseq+l[1]] - mat[l[0]*nseq+l[2]]) / 2;
    branch[2] = (mat[l[1]*nseq+l[2]] + mat[l[0]*nseq+l[2]] - mat[l[0]*nseq+l[1]]) / 2;

    branch[0] -= av[l[0]];
    branch[1] -= av[l[1]];
    branch[2] -= av[l[2]];

    // Print branch lengths
    //std::cout << branch[0] << std::endl;
    //std::cout << branch[1] << std::endl;
    //std::cout << branch[2] << std::endl;
    //std::cout << std::endl;

    // Reset negative branch lengths to zero
    if (branch[0] < 0) branch[0] = 0;
    if (branch[1] < 0) branch[1] = 0;
    if (branch[2] < 0) branch[2] = 0;

    // Build tree
    TVertexDescriptor internalVertex = seqan::addVertex(g);
    seqan::addEdge(g, internalVertex, seqan::getValue(connector, l[0]), (node_weight_type) roundToSignificantFigures((branch[0] / 10000000.0)* normFactor, 5));
    seqan::addEdge(g, internalVertex, seqan::getValue(connector, l[1]), (node_weight_type) roundToSignificantFigures((branch[1] / 10000000.0) * normFactor, 5));
    TVertexDescriptor the_root = seqan::addVertex(g);
    seqan::addEdge(g, the_root, seqan::getValue(connector, l[2]), (node_weight_type) roundToSignificantFigures((branch[2] / 20000000.0) * normFactor, 5));
    seqan::addEdge(g, the_root, internalVertex, (node_weight_type) roundToSignificantFigures((branch[2] / 20000000.0) * normFactor, 5));
    g.data_root = the_root;

    return g;
}

