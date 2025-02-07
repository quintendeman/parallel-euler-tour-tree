#pragma once

#include <utility>
#include <sequence/parallel_skip_list/include/skip_list_base.hpp>
#include <cassert>
#include <utilities/include/utils.h>


using std::pair;
constexpr int NA{-1};

namespace parallel_skip_list {

// Batch-parallel augmented skip list. Currently, the augmentation is
// hardcoded to the sum function with the value 1 assigned to each element. As
// such, `GetSum()` returns the size of the list.
//
// TODO(tomtseng): Allow user to pass in their own arbitrary associative
// augmentation functions. The contract for `GetSum` on a cyclic list should be
// that the function will be applied starting from `this`, because where we
// begin applying the function matters for non-commutative functions.
template <typename T>
class AugmentedElement : public ElementBase<AugmentedElement<T>> {
  friend class ElementBase<AugmentedElement<T>>;
 public:
  using ElementBase<AugmentedElement<T>>::Initialize;
  using ElementBase<AugmentedElement<T>>::Finish;

  static concurrent_array_allocator::Allocator<T>* val_allocator;

  static T* AllocateValueArray(int len) {
    T* values{val_allocator->Allocate(len)};
    for (int i = 0; i < len; i++) {
      values[i] = default_value;
    }
    return values;
  }

  static std::function<T(T,T)> aggregate_function;
  static T default_value;

  // See comments on `ElementBase<>`.
  AugmentedElement();
  explicit AugmentedElement(size_t random_int);
  ~AugmentedElement();

  // For each `{left, right}` in the `len`-length array `joins`, concatenate the
  // list that `left` lives in to the list that `right` lives in.
  //
  // `left` must be the last node in its list, and `right` must be the first
  // node of in its list. Each `left` must be unique, and each `right` must be
  // unique.
  static void BatchJoin(
      std::pair<AugmentedElement*, AugmentedElement*>* joins, int len);

  // For each `v` in the `len`-length array `splits`, split `v`'s list right
  // after `v`.
  static void BatchSplit(AugmentedElement** splits, int len);

  // For each `i`=0,1,...,`len`-1, assign value `new_values[i]` to element
  // `elements[i]`.
  static void BatchUpdate(AugmentedElement** elements, T* new_values, int len);

  // Assign value `new_value` to element `element`.
  static void Update(AugmentedElement* element, T new_value, int level = 0);

  // Get the result of applying the augmentation function over the subsequence
  // between `left` and `right` inclusive.
  //
  // `left` and `right` must live in the same list, and `left` must precede
  // `right` in the list.
  //
  // This function does not modify the data structure, so it may run
  // concurrently with other `GetSubsequenceSum` calls and const function calls.
  static int GetSubsequenceSum(
      const AugmentedElement* left, const AugmentedElement* right);

  // Get result of applying the augmentation function over the whole list that
  // the element lives in.
  T GetSum() const;

  using ElementBase<AugmentedElement>::FindRepresentative;
  using ElementBase<AugmentedElement>::GetPreviousElement;
  using ElementBase<AugmentedElement>::GetNextElement;

 private:
  static void DerivedInitialize();
  static void DerivedFinish();

  // Update aggregate value of node and clear `join_update_level` after joins.
  void UpdateTopDown(int level);
  void UpdateTopDownHelper(int level, AugmentedElement* curr);
  void UpdateTopDownSequential(int level);

  T* values_;
  // When updating augmented values, this marks the lowest index at which the
  // `values_` needs to be updated.
  int update_level_;
};

template<typename T>
concurrent_array_allocator::Allocator<T>* AugmentedElement<T>::val_allocator;

template<typename T>
std::function<T(T,T)> AugmentedElement<T>::aggregate_function;

template<typename T>
T AugmentedElement<T>::default_value;

template<typename T>
void AugmentedElement<T>::DerivedInitialize() {
  if (val_allocator == nullptr) {
    val_allocator = new concurrent_array_allocator::Allocator<T>;
  }
}

template<typename T>
void AugmentedElement<T>::DerivedFinish() {
  if (val_allocator != nullptr) {
    delete val_allocator;
  }
}

template<typename T>
AugmentedElement<T>::AugmentedElement() :
  ElementBase<AugmentedElement>{}, update_level_{NA} {
  values_ = AllocateValueArray(this->height_);
}

template<typename T>
AugmentedElement<T>::AugmentedElement(size_t random_int) :
  ElementBase<AugmentedElement>{random_int}, update_level_{NA} {
  values_ = AllocateValueArray(this->height_);
}

template<typename T>
AugmentedElement<T>::~AugmentedElement() {
  val_allocator->Free(values_, this->height_);
}

template<typename T>
void AugmentedElement<T>::UpdateTopDownSequential(int level) {
  if (level == 0) {
    if (this->height_ == 1) {
      update_level_ = NA;
    }
    return;
  }

  if (update_level_ < level) {
    UpdateTopDownSequential(level - 1);
  }
  int sum{values_[level - 1]};
  AugmentedElement* curr{this->neighbors_[level - 1].next};
  while (curr != nullptr && curr->height_ < level + 1) {
    if (curr->update_level_ != NA && curr->update_level_ < level) {
      curr->UpdateTopDownSequential(level - 1);
    }
    sum = aggregate_function(sum, curr->values_[level-1]);
    curr = curr->neighbors_[level - 1].next;
  }
  values_[level] = sum;

  if (this->height_ == level + 1) {
    update_level_ = NA;
  }
}

// `v.UpdateTopDown(level)` updates the augmented values of descendants of `v`'s
// `level`-th node. `update_level_` is used to determine what nodes need
// updating. `update_level_` is reset to `NA` for all traversed nodes at end of
// this function.
template<typename T>
void AugmentedElement<T>::UpdateTopDown(int level) {
  if (level <= 6) {
    UpdateTopDownSequential(level);
    return;
  }

  // Recursively update augmented values of children.
  UpdateTopDownHelper(level, this);

  // Now that children have correct augmented valeus, update self's augmented
  // value.
  int sum{values_[level - 1]};
  AugmentedElement* curr = this->neighbors_[level - 1].next;
  while (curr != nullptr && curr->height_ < level + 1) {
    sum = aggregate_function(sum, curr->values_[level-1]);
    curr = curr->neighbors_[level - 1].next;
  }
  values_[level] = sum;

  if (this->height_ == level + 1) {
    update_level_ = NA;
  }
}

template<typename T>
void AugmentedElement<T>::UpdateTopDownHelper(int level, AugmentedElement* curr) {
  if (curr->update_level_ != NA && curr->update_level_ < level) {
    parlay::parallel_do(
      [&] {
        auto next = curr->neighbors_[level-1].next;
        if (next != nullptr && next->height_ < level+1)
          UpdateTopDownHelper(level, next);
      },
      [&] { curr->UpdateTopDown(level-1); }
    );
  } else {
    auto next = curr->neighbors_[level-1].next;
    if (next != nullptr && next->height_ < level+1)
      UpdateTopDownHelper(level, next);
  }
}

// If `new_values` is non-null, for each `i`=0,1,...,`len`-1, assign value
// `new_vals[i]` to element `elements[i]`.
//
// If `new_values` is null, update the augmented values of the ancestors of
// `elements`, where the "ancestors" of element `v` refer to `v`,
// `v->FindLeftParent(0)`, `v->FindLeftParent(0)->FindLeftParent(1)`,
// `v->FindLeftParent(0)->FindLeftParent(2)`, and so on. This functionality is
// used privately to keep the augmented values correct when the list has
// structurally changed.
template<typename T>
void AugmentedElement<T>::BatchUpdate(
    AugmentedElement** elements, T* new_values, int len) {
  if (new_values != nullptr) {
    parallel_for (0, len, [&] (size_t i) {
      elements[i]->values_[0] = new_values[i];
    });
  }

  // The nodes whose augmented values need updating are the ancestors of
  // `elements`. Some nodes may share ancestors. `top_nodes` will contain,
  // without duplicates, the set of all ancestors of `elements` with no left
  // parents. From there we can walk down from those ancestors to update all
  // required augmented values.
  AugmentedElement** top_nodes{pbbs::new_array_no_init<AugmentedElement*>(len)};

  parallel_for (0, len, [&] (size_t i) {
    int level{0};
    AugmentedElement* curr{elements[i]};
    while (true) {
      int curr_update_level{curr->update_level_};
      if (curr_update_level == NA && CAS(&curr->update_level_, NA, level)) {
        level = curr->height_ - 1;
        AugmentedElement* parent{curr->FindLeftParent(level)};
        if (parent == nullptr) {
          top_nodes[i] = curr;
          break;
        } else {
          curr = parent;
          level++;
        }
      } else {
        // Someone other execution is shares this ancestor and has already
        // claimed it, so there's no need to walk further up.
        writeMin(&curr->update_level_, level);
        top_nodes[i] = nullptr;
        break;
      }
    }
  });

  parallel_for (0, len, [&] (size_t i) {
    if (top_nodes[i] != nullptr) {
      top_nodes[i]->UpdateTopDown(top_nodes[i]->height_ - 1);
    }
  });

  pbbs::delete_array(top_nodes, len);
}

template<typename T>
void AugmentedElement<T>::Update(AugmentedElement* element, T new_value, int level) {
  element->values_[level] = new_value;
  AugmentedElement* parent{element->FindLeftParent(level)};
  int sum{parent->values_[level]};
  AugmentedElement* curr{parent->neighbors_[level].next};
  while (curr != nullptr && curr->height_ == level+1) {
    sum = aggregate_function(sum, curr->values_[level]);
    curr = curr->neighbors_[level].next;
  }
  Update(parent, sum, level+1);
}

template<typename T>
void AugmentedElement<T>::BatchJoin(
    pair<AugmentedElement*, AugmentedElement*>* joins, int len) {
  AugmentedElement** join_lefts{
      pbbs::new_array_no_init<AugmentedElement*>(len)};
  parallel_for (0, len, [&] (size_t i) {
    AugmentedElement<T>::Join(joins[i].first, joins[i].second);
    join_lefts[i] = joins[i].first;
  });
  BatchUpdate(join_lefts, nullptr, len);
  pbbs::delete_array(join_lefts, len);
}

template<typename T>
void AugmentedElement<T>::BatchSplit(AugmentedElement** splits, int len) {
  parallel_for (0, len, [&] (size_t i) {
    splits[i]->Split();
  });
  parallel_for (0, len, [&] (size_t i) {
    AugmentedElement* curr{splits[i]};
    // `can_proceed` breaks ties when there are duplicate splits. When two
    // splits occur at the same place, only one of them should walk up and
    // update.
    bool can_proceed{
        curr->update_level_ == NA && CAS(&curr->update_level_, NA, 0)};
    if (can_proceed) {
      // Update values of `curr`'s ancestors.
      int sum{curr->values_[0]};
      int level{0};
      while (true) {
        if (level < curr->height_ - 1) {
          level++;
          curr->values_[level] = sum;
        } else {
          curr = curr->neighbors_[level].prev;
          if (curr == nullptr) {
            break;
          } else {
            sum += curr->values_[level];
          }
        }
      }
    }
  });
  parallel_for (0, len, [&] (size_t i) {
    splits[i]->update_level_ = NA;
  });
}

template<typename T>
int AugmentedElement<T>::GetSubsequenceSum(
    const AugmentedElement* left, const AugmentedElement* right) {
  int level{0};
  int sum{right->values_[level]};
  while (left != right) {
    level = min(left->height_, right->height_) - 1;
    if (level == left->height_ - 1) {
      sum += left->values_[level];
      left = left->neighbors_[level].next;
    } else {
      right = right->neighbors_[level].prev;
      sum += right->values_[level];
    }
  }
  return sum;
}

template<typename T>
T AugmentedElement<T>::GetSum() const {
  // Here we use knowledge of the implementation of `FindRepresentative()`.
  // `FindRepresentative()` gives some element that reaches the top level of the
  // list. For acyclic lists, the element is the leftmost one.
  AugmentedElement* root{FindRepresentative()};
  // Sum the values across the top level of the list.
  int level{root->height_ - 1};
  int sum{root->values_[level]};
  AugmentedElement* curr{root->neighbors_[level].next};
  while (curr != nullptr && curr != root) {
    sum += curr->values_[level];
    curr = curr->neighbors_[level].next;
  }
  if (curr == nullptr) {
    // The list is not circular, so we need to traverse backwards to beginning
    // of list and sum values to the left of `root`.
    curr = root;
    while (true) {
      while (level >= 0 && curr->neighbors_[level].prev == nullptr) {
        level--;
      }
      if (level < 0) {
        break;
      }
      while (curr->neighbors_[level].prev != nullptr) {
        curr = curr->neighbors_[level].prev;
        sum += curr->values_[level];
      }
    }
  }
  return sum;
}

}  // namespace parallel_skip_list
