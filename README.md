# ParlayLib Version of Batch-Parallel Euler Tour Trees

This code is forked from the code for the ALENEX 2019 paper [Batch-Parallel Euler Tour Trees](https://arxiv.org/abs/1810.10738).
The original README.md file is appended below for.

In this repository I converted the code to use the parallel functionality of the ParlayLib library.
I only converted the parallel skip list and parallel Euler tour tree code, not the parallel treaps code or any sequential code.
I also added support for a custom subtree value type and augmentation, which needs to be tested.

## Including the ParlayLib version through CMake

This excerpt shows an example of how to include the Batch-Parallel Euler Tour Trees code in this repository using CMake.

```
FetchContent_Declare(parallel-ett
    GIT_REPOSITORY  git@github.com:quintendeman/parallel-euler-tour-tree.git
    GIT_TAG         master
)
FetchContent_MakeAvailable(parallel-ett)

add_executable(test
  test/test.cpp
)
target_link_libraries(test PRIVATE par-ett)
```

## Using the Library

This small excerpt shows an example of using the library once included through CMake.

```
#include "ParETT/euler_tour_tree.hpp"

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
```

## Using Custom Augmentation

The code will still use integer values and the sum function by default.
To customize the augmentations use the desired type in the template.
Then define the two static variables similar to the example here:

```
parallel_skip_list::AugmentedElement<int>::aggregate_function = [&] (int x, int y) { return x + y; };
parallel_skip_list::AugmentedElement<int>::default_value = 1;
```


# Batch-parallel Euler tour trees

This is the code for the ALENEX 2019 paper [Batch-Parallel Euler Tour
Trees](https://arxiv.org/abs/1810.10738).

The items of interest in this repository are the implementations of the
batch-parallel skip list (`src/sequence/parallel_skip_list`) and the
batch-parallel Euler tour tree (`src/dynamic_trees/parallel_euler_tour_tree`).
Other folders in `src/sequence/` and `src/dynamic_trees/` contain code used for
benchmarks in the paper.

## Compilation

This code is written assuming that the compiler is g++ 5.5.0 with Cilk Plus
extensions.

Existing benchmarking and testing code may be compiled by using the Makefiles
scattered throughout `src/`. Executables resulting from compilation are placed
in `bin/`.

## Data structure descriptions

### Skip lists

In `src/sequence/parallel_skip_list`, we implement both augmented and
unaugmented batch-parallel skip lists supporting _joining_ (concatenating two
sequences) and _splitting_ (cleaving a sequence into two at a chosen point).
The skip lists can execute a batch of _k_ joins, _k_ splits, or _k_ augmented
value updates over _n_ elements in _O(k log(1 + n/k))_ expected work and _O(log
n)_ depth with high probability.

In addition, the unaugmented skip lists are _phase-concurrent_. Informally, what
this means that instead of having to organize calls into a batch and feed them
all in at once into `BatchJoin` (or `BatchSplit`), each thread can individually
call `Join` (or `Split`) asynchronously.

### Euler tour trees

The Euler tour tree is a data structure for the dynamic trees problem. In the
dynamic trees problem, we want to maintain a forest undergoing _links_ (edge
additions) and _cuts_ (edge deletions). We also want to be able to query whether
pairs of vertices are connected.

In our paper, we describe a batch-parallel Euler tour tree that can execute a
batch of _k_ links or _k_ cuts over an _n_-vertex forest in _O(k
log(1 + n/k))_ expected work and _O(log n)_ depth with high probability. We
implement this with some simplifying modifications in
`src/dynamic_trees/parallel_euler_tour_tree`.

## Future work on this repository
* Currently the augmented skip list only augments with the size of the list.
  We should allow the user to specify the augmentation function.
* We should try building Euler tour trees on top of augmented skip lists instead
  of unaugmented ones and support querying of sizes of trees in the represented
  forest.
* There are at most _3n - 2_ elements in an Euler tour tree at any given time. Can
  we get noticeably better performance by preallocating these elements at
  initialization and reusing them instead of repeatedly allocating and
  deallocating them?
* Better tests should be written, and the tests should be migrated to Google
  Test or another nice testing framework.
* The parallel skip list and parallel Euler tour tree code has been cleaned up,
  but all the rest of the code in this repository is in a poorer state.
* The build system consists of lots of Makefiles. The Makefiles perform
  in-source builds and put executables in bin/ (for easy gitignoring) with no
  subdirectories.  It would be nice to have out-of-source builds and to give
  `bin/` a directory structure matching that of `src/` so that we don't have to
  worry about different executables with the same name overwriting each other.
  Switching to CMake would make these changes easier.
* Proper `make clean` commands are not implemented for the benchmarks.

## Resources
Thomas Tseng, Laxman Dhulipala, and Guy Blelloch. [Batch-parallel Euler tour
trees](https://arxiv.org/abs/1810.10738). In _Proceedings of the Twenty-First
Workshop on Algorithm Engineering and Experiments_, page to appear.  Society for
Industrial and Applied Mathematics, 2019.
