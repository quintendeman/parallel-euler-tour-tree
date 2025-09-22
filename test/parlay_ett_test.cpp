#include <gtest/gtest.h>
#include "dynamic_trees/parallel_euler_tour_tree/include/euler_tour_tree.hpp"
#include "dynamic_trees/parallel_euler_tour_tree/include/unaugmented_euler_tour_tree.hpp"


TEST(ParlaySuite, mini_unaugmented_test) {
    int n = 1000;
    int k = n-1;
    int num_trials = 100;
    srand(time(NULL));

    using EulerTourTree = parallel_euler_tour_tree::UnaugmentedEulerTourTree;

    for (int trial = 0; trial < num_trials; trial++) {
        size_t seed = rand();
        EulerTourTree tree(n, seed);
        parlay::sequence<std::pair<int,int>> links(k);
        for (int i = 0; i < k; i++)
            links[i] = {i,i+1};
        tree.BatchLink(links);
        tree.BatchCut(links);
    }
}

TEST(ParlaySuite, mini_aggregate_test) {
    int n = 1000;
    int k = n-1;
    int num_trials = 100;
    srand(time(NULL));

    using EulerTourTree = parallel_euler_tour_tree::EulerTourTree<int>;
    parallel_skip_list::AugmentedElement<int>::default_value = 1;
    parallel_skip_list::AugmentedElement<int>::aggregate_function = [] (int x, int y) { return x+y; };

    for (int trial = 0; trial < num_trials; trial++) {
        size_t seed = rand();
        EulerTourTree tree(n, seed);
        parlay::sequence<std::pair<int,int>> links(k);
        for (int i = 0; i < k; i++)
            links[i] = {i,i+1};
        tree.BatchLink(links);
        ASSERT_EQ(tree.vertices_[0].GetSum(), n + 2*k) << "INCORRECT AGGREGATE AFTER BATCH LINK." << std::endl;
        tree.BatchCut(links);
        ASSERT_EQ(tree.vertices_[0].GetSum(), 1) << "INCORRECT AGGREGATE AFTER BATCH CUT." << std::endl;
    }
}
