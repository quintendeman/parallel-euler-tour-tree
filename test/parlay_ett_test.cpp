#include <gtest/gtest.h>
#include "dynamic_trees/parallel_euler_tour_tree/include/euler_tour_tree.hpp"


TEST(ParlaySuite, test) {
    int n = 1000;
    int k = 250;
    using EulerTourTree = parallel_euler_tour_tree::EulerTourTree<int>;
    EulerTourTree tree(n);
    std::pair<int,int>* links = new std::pair<int,int>[k];
    for (int i = 0; i < k; i++)
        links[i] = {i,i+1};
    tree.BatchLink(links, k);
    tree.BatchCut(links, k);
    delete[] links;
}
