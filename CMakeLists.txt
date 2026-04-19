cmake_minimum_required(VERSION 3.20)
project(cancel_tree_majority LANGUAGES CXX)

enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(MOCKTURTLE_ROOT "/home/mrunal/Desktop/mockturtle" CACHE PATH "Path to mockturtle root")

if(NOT EXISTS "${MOCKTURTLE_ROOT}/include/mockturtle")
  message(FATAL_ERROR "Could not find mockturtle headers under ${MOCKTURTLE_ROOT}/include/mockturtle")
endif()

add_library(majlib INTERFACE)

target_include_directories(majlib INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${MOCKTURTLE_ROOT}/include
  ${MOCKTURTLE_ROOT}/lib/fmt
  ${MOCKTURTLE_ROOT}/lib/kitty
  ${MOCKTURTLE_ROOT}/lib/lorina
  ${MOCKTURTLE_ROOT}/lib/parallel_hashmap
  ${MOCKTURTLE_ROOT}/lib/rang
  ${MOCKTURTLE_ROOT}/lib/json
)

target_compile_options(majlib INTERFACE
  -Wall
  -Wextra
  -Wpedantic
)

target_compile_definitions(majlib INTERFACE
  FMT_HEADER_ONLY=1
)

add_executable(majgen src/majgen.cpp)
target_link_libraries(majgen PRIVATE majlib)

add_executable(test_cancel_tree test/test_cancel_tree.cpp)
target_link_libraries(test_cancel_tree PRIVATE majlib)

add_test(NAME cancel_tree_correctness COMMAND test_cancel_tree)
