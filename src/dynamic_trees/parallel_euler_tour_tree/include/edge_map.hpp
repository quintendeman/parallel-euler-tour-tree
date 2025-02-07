#pragma once

#include <utility>
#include <parlay/parallel.h>
#include <parlay/alloc.h>

#include <utilities/include/concurrentMap.h>
#include <utilities/include/hash_pair.hpp>
#include <dynamic_trees/parallel_euler_tour_tree/include/euler_tour_sequence.hpp>

namespace parallel_euler_tour_tree {

namespace _internal {

// Used in Euler tour tree for mapping directed edges (pairs of ints) to the
// sequence element in the Euler tour representing the edge.
//
// Only one of (u, v) and (v, u) should be added to the map; we can find the
// other edge using the `twin_` pointer in `Element`.
template<typename T>
class EdgeMap {
using Element = _internal::Element<T>;
 public:
  EdgeMap() = delete;
  explicit EdgeMap(int num_vertices);
  ~EdgeMap();

  bool Insert(int u, int v, Element* edge);
  bool Delete(int u, int v);
  Element* Find(int u, int v);

  // Deallocate all elements held in the map. This assumes that all elements
  // in the map were allocated through `allocator`.
  void FreeElements(parlay::type_allocator<Element>* allocator);

 private:
  concurrent_map::concurrentHT<
      std::pair<int, int>, Element*, HashIntPairStruct> map_;
};


template<typename T>
EdgeMap<T>::EdgeMap(int num_vertices)
    : map_{nullptr, static_cast<size_t>(num_vertices - 1),
          std::make_pair(-1, -1), std::make_pair(-2, -2)} {}

template<typename T>
EdgeMap<T>::~EdgeMap() {
  map_.del();
}

template<typename T>
bool EdgeMap<T>::Insert(int u, int v, Element* edge) {
  if (u > v) {
    std::swap(u, v);
    edge = edge->twin_;
  }
  return map_.insert(make_pair(u, v), edge);
}

template<typename T>
bool EdgeMap<T>::Delete(int u, int v) {
  if (u > v) {
    std::swap(u, v);
  }
  return map_.deleteVal(make_pair(u, v));
}

template<typename T>
_internal::Element<T>* EdgeMap<T>::Find(int u, int v) {
  if (u > v) {
    Element* vu{*map_.find(make_pair(v, u))};
    return vu == nullptr ? nullptr : vu->twin_;
  } else {
    return *map_.find(make_pair(u, v));
  }
}

template<typename T>
void EdgeMap<T>::FreeElements(parlay::type_allocator<Element>* allocator) {
  parallel_for (0, map_.capacity, [&] (size_t i) {
    auto kv{map_.table[i]};
    auto key{get<0>(kv)};
    if (key != map_.empty_key && key != map_.tombstone) {
      Element* element{get<1>(kv)};
      element->twin_->~Element();
      allocator->free(element->twin_);
      element->~Element();
      allocator->free(element);
    }
  });
}

}  // namespace _internal

}  // namespace parallel_euler_tour_tree
