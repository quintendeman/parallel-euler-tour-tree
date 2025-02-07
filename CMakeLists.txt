cmake_minimum_required(VERSION 3.15)
project(ParETT)
include (FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)
# add_compile_options(-fsanitize=undefined)
# add_link_options(-fsanitize=undefined)

FetchContent_Declare(parlaylib
  GIT_REPOSITORY  https://github.com/cmuparlay/parlaylib.git
  GIT_TAG         master
)
FetchContent_GetProperties(parlaylib)
if(NOT parlaylib_POPULATED)
  FetchContent_Populate(parlaylib)  
  add_subdirectory(${parlaylib_SOURCE_DIR} EXCLUDE_FROM_ALL)
endif()

# FetchContent_Declare(abseil-cpp
#   GIT_REPOSITORY  https://github.com/abseil/abseil-cpp.git
#   GIT_TAG         master
# )
# FetchContent_MakeAvailable(abseil-cpp)

add_library(par-ett INTERFACE)
include(GNUInstallDirs)
set(PARETT_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(par-ett INTERFACE
  $<BUILD_INTERFACE:${PARETT_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  src
)

add_executable(skip_list_tests
  src/sequence/parallel_skip_list/tests/test_parallel_skip_list.cpp
)
target_link_libraries(skip_list_tests PRIVATE parlay)
target_include_directories(skip_list_tests PRIVATE src)

add_executable(augmented_skip_list_tests
  src/sequence/parallel_skip_list/tests/test_parallel_augmented_skip_list.cpp
)
target_link_libraries(augmented_skip_list_tests PRIVATE parlay)
target_include_directories(augmented_skip_list_tests PRIVATE src)

add_executable(euler_tour_tree_tests
  src/dynamic_trees/parallel_euler_tour_tree/tests/test_parallel_euler_tour_tree.cpp
  src/dynamic_trees/parallel_euler_tour_tree/tests/simple_forest_connectivity.cpp
)
target_link_libraries(euler_tour_tree_tests PRIVATE parlay)
target_include_directories(euler_tour_tree_tests PRIVATE src)

add_executable(benchmark_skip_list
  src/sequence/benchmarks/batch_sequence_benchmark/parallel_augmented_skip_list/benchmark_batch_sequence_parallel_augmented_skip_list.cpp
)
target_link_libraries(benchmark_skip_list PRIVATE parlay)
target_include_directories(benchmark_skip_list PRIVATE src)

add_executable(benchmark_euler_tour_tree
  src/dynamic_trees/benchmarks/parallel_ett/benchmark_dynamic_trees_parallel_ett.cpp
)
target_link_libraries(benchmark_euler_tour_tree PRIVATE parlay)
target_include_directories(benchmark_euler_tour_tree PRIVATE src)