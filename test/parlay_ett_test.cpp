#include <gtest/gtest.h>
#include "dynamic_trees/parallel_euler_tour_tree/include/euler_tour_tree.hpp"


TEST(ParlaySuite, test) {
    int num_trials = 10;
    int n = 1000;
    int k = 250;
    using EulerTourTree = parallel_euler_tour_tree::EulerTourTree<int>;

    for (int trial = 0; trial < num_trials; trial++) {
        EulerTourTree tree(n);
        parlay::sequence<std::pair<int,int>> links(k);
        for (int i = 0; i < k; i++)
            links[i] = {i,i+1};
        tree.BatchLink(links);
        tree.BatchCut(links);
    }
}
