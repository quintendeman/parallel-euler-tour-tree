#pragma once

#include <utility>

#include <dynamic_trees/parallel_euler_tour_tree/include/edge_map.hpp>
#include <dynamic_trees/parallel_euler_tour_tree/include/euler_tour_sequence.hpp>
#include <sequence/parallel_skip_list/include/skip_list_base.hpp>
#include <utilities/include/random.h>
#include <utilities/include/seq.h>
#include <utilities/include/sequence_ops.h>
#include <utilities/include/utils.h>

namespace parallel_euler_tour_tree {

// Euler tour trees represent forests. We may add an edge using `Link`, remove
// an edge using `Cut`, and query whether two vertices are in the same tree
// using `IsConnected`. This implementation can also exploit parallelism when
// many edges are added at once through `BatchLink` or many edges are deleted at
// once through `BatchCut`.
template<typename T>
class EulerTourTree {
using Element = _internal::Element<T>;
 public:
  static parlay::type_allocator<_internal::Element<T>> allocator;

  EulerTourTree() = delete;
  // Initializes n-vertex forest with no edges.
  explicit EulerTourTree(int num_vertices);
  ~EulerTourTree();
  EulerTourTree(const EulerTourTree&) = delete;
  EulerTourTree(EulerTourTree&&) = delete;
  EulerTourTree& operator=(const EulerTourTree&) = delete;
  EulerTourTree& operator=(EulerTourTree&&) = delete;

  // Returns true if `u` and `v` are in the same tree in the represented forest.
  bool IsConnected(int u, int v) const;
  // Uses `parent` pointers to find representatives faster.
  bool IsConnected2(int u, int v) const;
  // Adds edge {`u`, `v`} to forest. The addition of this edge must not create a
  // cycle in the graph.
  void Link(int u, int v);
  void Link2(int u, int v);
  void Link3(int u, int v);
  // Removes edge {`u`, `v`} from forest. The edge must be present in the
  // forest.
  void Cut(int u, int v);
  // Set the value of vertex `v` to be `new_value`.
  void Update(int v, T new_value);
  // Update the augmented value of vertex `v` and each node that aggregates it by
  // applying the function `f` on the value of each node.
  void UpdateWithFunction(int v, std::function<void(T&)> f);
  // Return the value associated with element `v`.
  T GetValue(int v);

  // Adds all edges in the `len`-length array `links` to the forest. Adding
  // these edges must not create cycles in the graph.
  void BatchLink(std::pair<int, int>* links, int len);
  // Removes all edges in the `len`-length array `cuts` from the forest. These
  // edges must be present in the forest and must be distinct.
  void BatchCut(std::pair<int, int>* cuts, int len);
  // Updates all the vertices in the `len`-length array `vertices` with the
  // new corresponding value in the `new_values` array.
  void BatchUpdate(int* vertices, T* new_values, int len);

 private:
  void BatchCutRecurse(std::pair<int, int>* cuts, int len,
      bool* ignored, _internal::Element<T>** join_targets,
      _internal::Element<T>** edge_elements);

  int num_vertices_;
  pbbs::random randomness_;

  std::vector<Element*> node_pool;
 public:
  _internal::Element<T>* vertices_;
  _internal::EdgeMap<T> edges_;
};

template<typename T>
parlay::type_allocator<_internal::Element<T>> EulerTourTree<T>::allocator;



namespace {

  // On BatchCut, randomly ignore 1/`kBatchCutRecursiveFactor` cuts and recurse
  // on them later.
  constexpr int kBatchCutRecursiveFactor{100};

  template<typename T>
  void BatchCutSequential(EulerTourTree<T>* ett, pair<int, int>* cuts, int len) {
    for (int i = 0; i < len; i++) {
      ett->Cut(cuts[i].first, cuts[i].second);
    }
  }

  template<typename T>
  void BatchLinkSequential(EulerTourTree<T>* ett, pair<int, int>* links, int len) {
    for (int i = 0; i < len; i++) {
      ett->Link(links[i].first, links[i].second);
    }
  }

}  // namespace

template<typename T>
EulerTourTree<T>::EulerTourTree(int num_vertices)
    : num_vertices_{num_vertices} , edges_{num_vertices_} , randomness_{} {
  Element::Initialize();
  vertices_ = pbbs::new_array_no_init<Element>(num_vertices_);
  parallel_for (0, num_vertices_, [&] (size_t i) {
    new (&vertices_[i]) Element{randomness_.ith_rand(i)};
    // The Euler tour on a vertex v (a singleton tree) is simply (v, v).
    Element::Join(&vertices_[i], &vertices_[i]);
  });
  randomness_ = randomness_.next();
  for (int i = 0; i < 3*num_vertices-2; i++)
    node_pool.push_back(allocator.create(randomness_.ith_rand(i)));
}

template<typename T>
EulerTourTree<T>::~EulerTourTree() {
  pbbs::delete_array(vertices_, num_vertices_);
  for (auto node : node_pool)
    allocator.destroy(node);
  edges_.FreeElements(&allocator);
  Element::Finish();
}

template<typename T>
bool EulerTourTree<T>::IsConnected(int u, int v) const {
  return vertices_[u].FindRepresentative() == vertices_[v].FindRepresentative();
}

template<typename T>
bool EulerTourTree<T>::IsConnected2(int u, int v) const {
  return vertices_[u].FindRepresentative2() == vertices_[v].FindRepresentative2();
}

template<typename T>
void EulerTourTree<T>::Link(int u, int v) {
  Element* uv = node_pool.back();
  node_pool.pop_back();
  Element* vu = node_pool.back();
  node_pool.pop_back();
  randomness_ = randomness_.next();
  uv->twin_ = vu;
  vu->twin_ = uv;
  edges_.Insert(u, v, uv);
  Element* u_left{&vertices_[u]};
  Element* v_left{&vertices_[v]};
  Element* u_right = (Element*) u_left->SequentialSplitRight(false);
  Element* v_right = (Element*) v_left->SequentialSplitRight(false);
  Element::SequentialJoin(u_left, uv, false);
  Element::SequentialJoin(uv, v_right, false);
  Element::SequentialJoin(v_left, vu, false);
  Element::SequentialJoin(vu, u_right, false);
  Element::Update(u_left, u_left->values_[0]);
  Element::Update(v_left, v_left->values_[0]);
  Element::Update(u_right, u_right->values_[0]);
  Element::Update(v_right, v_right->values_[0]);
}

template<typename T>
void EulerTourTree<T>::Link2(int u, int v) {
  Element* uv = node_pool.back();
  node_pool.pop_back();
  Element* vu = node_pool.back();
  node_pool.pop_back();
  randomness_ = randomness_.next();
  uv->twin_ = vu;
  vu->twin_ = uv;
  edges_.Insert(u, v, uv);
  Element* u_right{&vertices_[u]};
  Element* v_right{&vertices_[v]};
  Element* u_left = (Element*) u_right->SequentialSplitLeft(false);
  Element* v_left = (Element*) v_right->SequentialSplitLeft(false);
  Element::SequentialJoin2(u_left, uv, false);
  Element::SequentialJoin2(uv, v_right, false);
  Element::SequentialJoin2(v_left, vu, false);
  Element::SequentialJoin2(vu, u_right, false);
  Element::Update(u_left, u_left->values_[0]);
  Element::Update(v_left, v_left->values_[0]);
  Element::Update(u_right, u_right->values_[0]);
  Element::Update(v_right, v_right->values_[0]);
}

template<typename T>
void EulerTourTree<T>::Link3(int u, int v) {
  Element* uv = node_pool.back();
  node_pool.pop_back();
  Element* vu = node_pool.back();
  node_pool.pop_back();
  randomness_ = randomness_.next();
  uv->twin_ = vu;
  vu->twin_ = uv;
  edges_.Insert(u, v, uv);
  Element* u_left{&vertices_[u]};
  Element* v_right{&vertices_[v]};
  Element* u_right = (Element*) u_left->SequentialSplitRight(false);
  Element* v_left = (Element*) v_right->SequentialSplitLeft(false);
  Element::SequentialJoin2(u_left, uv, false);
  Element::SequentialJoin2(uv, v_left, false);
  Element::SequentialJoin2(v_right, vu, false);
  Element::SequentialJoin2(vu, u_right, false);
  Element::Update(u_left, u_left->values_[0]);
  Element::Update(v_left, v_left->values_[0]);
  Element::Update(u_right, u_right->values_[0]);
  Element::Update(v_right, v_right->values_[0]);
}

template<typename T>
void EulerTourTree<T>::BatchLink(pair<int, int>* links, int len) {
  if (len <= 75) {
    BatchLinkSequential(this, links, len);
    return;
  }

  // For each added edge {x, y}, allocate elements (x, y) and (y, x).
  // For each vertex x that shows up in an added edge, split on (x, x). Let
  // succ(x) denote the successor of (x, x) prior to splitting.
  // For each vertex x, identify which y_1, y_2, ... y_k that x will be newly
  // connected to by performing a semisort on {(x, y), (y, x) : {x, y} is an
  // added edge}.
  // If x has new neighbors y_1, y_2, ..., y_k, join (x, x) to (x, y_1). Join
  // (y_i,x) to (x, y_{i+1}) for each i < k. Join (y_k, x) to succ(x).

  parlay::sequence<pair<uint32_t,uint32_t>> links_both_dirs(2*len);
  parallel_for (0, len, [&] (size_t i) {
    links_both_dirs[2 * i] = links[i];
    links_both_dirs[2 * i + 1] = make_pair(links[i].second, links[i].first);
  });
  parlay::integer_sort_inplace(links_both_dirs, [&] (pair<uint32_t,uint32_t> p) { return p.first; });

  Element** split_successors{pbbs::new_array_no_init<Element*>(2 * len)};
  parallel_for (0, 2*len, [&] (size_t i) {
    int u, v;
    std::tie(u, v) = links_both_dirs[i];

    // split on each vertex that appears in the input
    if (i == 2 * len - 1 || u != links_both_dirs[i + 1].first) {
      split_successors[i] = (Element*) vertices_[u].Split();
    }

    // allocate edge element
    if (u < v) {
      Element* uv{allocator.create(randomness_.ith_rand(2*i))};
      Element* vu{allocator.create(randomness_.ith_rand(2*i+1))};
      uv->twin_ = vu;
      vu->twin_ = uv;
      edges_.Insert(u, v, uv);
    }
  });
  randomness_ = randomness_.next();

  parallel_for (0, 2*len, [&] (size_t i) {
    int u, v;
    std::tie(u, v) = links_both_dirs[i];
    Element* uv{edges_.Find(u, v)};
    Element* vu{uv->twin_};
    if (i == 0 ||
        u != links_both_dirs[i - 1].first) {
      Element::Join(&vertices_[u], uv);
    }
    if (i == 2 * len - 1 ||
        u != links_both_dirs[i + 1].first) {
      Element::Join(vu, split_successors[i]);
    } else {
      int u2, v2;
      std::tie(u2, v2) = links_both_dirs[i + 1];
      Element::Join(vu, edges_.Find(u2, v2));
    }
  });

  pbbs::delete_array(split_successors, 2 * len);
}

template<typename T>
void EulerTourTree<T>::Cut(int u, int v) {
  Element* uv{edges_.Find(u, v)};
  Element* vu{uv->twin_};
  edges_.Delete(u, v);
  Element* u_left = (Element*) uv->GetPreviousElement();
  Element* v_left = (Element*) vu->GetPreviousElement();
  Element* v_right = (Element*) uv->SequentialSplitRight(false);
  Element* u_right = (Element*) vu->SequentialSplitRight(false);
  u_left->SequentialSplitRight(false);
  v_left->SequentialSplitRight(false);
  Element::Update(uv, uv->values_[0]);
  Element::Update(vu, vu->values_[0]);
  node_pool.push_back(uv);
  node_pool.push_back(vu);
  Element::SequentialJoin(u_left, u_right, false);
  Element::SequentialJoin(v_left, v_right, false);
  Element::Update(u_left, u_left->values_[0]);
  Element::Update(v_left, v_left->values_[0]);
  Element::Update(u_right, u_right->values_[0]);
  Element::Update(v_right, v_right->values_[0]);
}

// `ignored`, `join_targets`, and `edge_elements` are scratch space.
// `ignored[i]` will be set to true if `cuts[i]` will not be executed in this
// round of recursion.
// `join_targets` stores sequence elements that need to be joined to each other.
// `edge_elements[i]` stores a pointer to the sequence element corresponding to
// edge `cuts[i]`.
template<typename T>
void EulerTourTree<T>::BatchCutRecurse(pair<int, int>* cuts, int len,
    bool* ignored, Element** join_targets, Element** edge_elements) {
  if (len <= 75) {
    BatchCutSequential(this, cuts, len);
    return;
  }

  // Notation: "(x, y).next" is the next element in the tour (x, y) is in. "(x,
  // y).prev" is the previous element. "(x, y).twin" is (y, x).
  // For each edge {x, y} to cut:
  // Sequentially, we'd want to join (y, x).prev to (x, y).next and (x, y).prev
  // to (y, x).next. We can't correctly do this if any of those four elements
  // are to be cut and removed as well. Instead, for dealing with connecting (y,
  // x).prev to (x, y).next (dealing with connecting (x, y).prev to (y,x ).next
  // is symmetric), we do the following:
  // - If (y, x).prev is to be cut, then do nothing --- some other thread will
  // deal with this.
  // - Otherwise, start with element e = (x, y).next. So long as e is to be cut,
  // traverse to the next possible join location at e := e.next.twin. Join
  // (y,x).prev to e.
  // This strategy doesn't have good depth since we may have to traverse on e
  // for a long time. To fix this, we randomly ignore some cuts so that all
  // traversal lengths are O(log n) with high probability. We perform all
  // unignored cuts as described above, and recurse on the ignored cuts
  // afterwards.

  parallel_for (0, len, [&] (size_t i) {
    ignored[i] = randomness_.ith_rand(i) % kBatchCutRecursiveFactor == 0;

    if (!ignored[i]) {
      int u, v;
      std::tie(u, v) = cuts[i];
      Element* uv{edges_.Find(u, v)};
      edge_elements[i] = uv;
      Element* vu{uv->twin_};
      uv->split_mark_ = vu->split_mark_ = true;
    }
  });
  randomness_ = randomness_.next();

  parallel_for (0, len, [&] (size_t i) {
    if (!ignored[i]) {
      Element* uv{edge_elements[i]};
      Element* vu{uv->twin_};

      Element* left_target = (Element*) uv->GetPreviousElement();
      if (left_target->split_mark_) {
        join_targets[4 * i] = nullptr;
      } else {
        Element* right_target = (Element*) vu->GetNextElement();
        while (right_target->split_mark_) {
          right_target = (Element*) right_target->twin_->GetNextElement();
        }
        join_targets[4 * i] = left_target;
        join_targets[4 * i + 1] = right_target;
      }

      left_target = (Element*) vu->GetPreviousElement();
      if (left_target->split_mark_) {
        join_targets[4 * i + 2] = nullptr;
      } else {
        Element* right_target = (Element*) uv->GetNextElement();
        while (right_target->split_mark_) {
          right_target = (Element*) right_target->twin_->GetNextElement();
        }
        join_targets[4 * i + 2] = left_target;
        join_targets[4 * i + 3] = right_target;
      }
    }
  });

  parallel_for (0, len, [&] (size_t i) {
    if (!ignored[i]) {
      Element* uv{edge_elements[i]};
      Element* vu{uv->twin_};
      uv->Split();
      vu->Split();
      Element* predecessor = (Element*) uv->GetPreviousElement();
      if (predecessor != nullptr) {
        predecessor->Split();
      }
      predecessor = (Element*) vu->GetPreviousElement();
      if (predecessor != nullptr) {
        predecessor->Split();
      }
    }
  });

  parallel_for (0, len, [&] (size_t i) {
    if (!ignored[i]) {
      // Here we must use `edge_elements[i]` instead of `edges_.Find(u, v)`
      // because the concurrent hash table cannot handle simultaneous lookups
      // and deletions.
      Element* uv{edge_elements[i]};
      Element* vu{uv->twin_};
      allocator.destroy(uv);
      allocator.destroy(vu);
      int u, v;
      std::tie(u, v) = cuts[i];
      edges_.Delete(u, v);

      if (join_targets[4 * i] != nullptr) {
        Element::Join(join_targets[4 * i], join_targets[4 * i + 1]);
      }
      if (join_targets[4 * i + 2] != nullptr) {
        Element::Join(join_targets[4 * i + 2], join_targets[4 * i + 3]);
      }
    }
  });

  seq::sequence<pair<int, int>> cuts_seq{
      seq::sequence<pair<int, int>>(cuts, len)};
  seq::sequence<bool> ignored_seq{seq::sequence<bool>(ignored, len)};
  seq::sequence<pair<int, int>> next_cuts_seq{
    pbbs::pack(cuts_seq, ignored_seq)};
  BatchCutRecurse(next_cuts_seq.as_array(), next_cuts_seq.size(),
      ignored, join_targets, edge_elements);
  pbbs::delete_array(next_cuts_seq.as_array(), next_cuts_seq.size());
}

template<typename T>
void EulerTourTree<T>::BatchCut(pair<int, int>* cuts, int len) {
  if (len <= 75) {
    BatchCutSequential(this, cuts, len);
    return;
  }
  bool* ignored{pbbs::new_array_no_init<bool>(len)};
  Element** join_targets{pbbs::new_array_no_init<Element*>(4 * len)};
  Element** edge_elements{pbbs::new_array_no_init<Element*>(len)};
  BatchCutRecurse(cuts, len, ignored, join_targets, edge_elements);
  pbbs::delete_array(edge_elements, len);
  pbbs::delete_array(join_targets, 4 * len);
  pbbs::delete_array(ignored, len);
}

template<typename T>
void EulerTourTree<T>::Update(int v, T new_value) {
  Element::Update(&vertices_[v], new_value);
}

template<typename T>
void EulerTourTree<T>::UpdateWithFunction(int v, std::function<void(T&)> f) {
  Element::UpdateWithFunction(&vertices_[v], f);
}

template<typename T>
void EulerTourTree<T>::BatchUpdate(int* vertices, T* new_values, int len) {
  Element** update_targets{pbbs::new_array_no_init<Element*>(len)};
  parallel_for (0, len, [&] (size_t i) {
    update_targets[i] = vertices_[vertices[i]];
  });
  BatchUpdate(update_targets, new_values, len);
}

template<typename T>
T EulerTourTree<T>::GetValue(int v) {
  return vertices_[v].values_[0];
}

}  // namespace parallel_euler_tour_tree
