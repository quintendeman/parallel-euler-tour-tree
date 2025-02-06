#include <sequence/parallel_skip_list/include/augmented_skip_list.hpp>

#include <string>

#include <sequence/benchmarks/batch_sequence_benchmark/benchmark.hpp>
#include <utilities/include/utils.h>

namespace bsb = batch_sequence_benchmark;
using Element = parallel_skip_list::AugmentedElement<int>;
using std::string;

int main(int argc, char** argv) {
  Element::aggregate_function = [&] (int x, int y) { return x + y; };
  Element::default_value = 1;
  bsb::BenchmarkParameters parameters{bsb::GetBenchmarkParameters(argc, argv)};
  Element::Initialize();
  Element* elements{pbbs::new_array_no_init<Element>(parameters.num_elements)};
  pbbs::random r{};
  parlay::parallel_for (0, parameters.num_elements, [&] (size_t i) {
    new (&elements[i]) Element{r.ith_rand(i)};
  });

  bsb::RunBenchmark(elements, parameters);

  pbbs::delete_array(elements, parameters.num_elements);
  Element::Finish();
}
