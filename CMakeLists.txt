cmake_minimum_required(VERSION 3.15)
project(ParallelETT)
include (FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -march=native")

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

add_executable(skip_list_tests
  src/sequence/parallel_skip_list/tests/test_parallel_skip_list.cpp
  src/sequence/parallel_skip_list/src/skip_list_base.cpp
)
target_link_libraries(skip_list_tests PRIVATE parlay)
target_include_directories(skip_list_tests PRIVATE src)