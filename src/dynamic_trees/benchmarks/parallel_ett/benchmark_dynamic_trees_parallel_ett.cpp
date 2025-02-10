#include <dynamic_trees/parallel_euler_tour_tree/include/euler_tour_tree.hpp>

#include <dynamic_trees/benchmarks/benchmark.hpp>

int main(int argc, char** argv) {
  parallel_skip_list::AugmentedElement<int>::aggregate_function = [&] (int x, int y) { return x + y; };
  parallel_skip_list::AugmentedElement<int>::default_value = 1;
  dynamic_trees_benchmark::RunBenchmark<
      parallel_euler_tour_tree::EulerTourTree<int>>(argc, argv);
  return 0;
}
