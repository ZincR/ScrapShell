cmake_minimum_required(VERSION 3.20)
project(phcpf_mesher CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(DIRECTIONAL_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../third_party/Directional)

# GMP is required: Directional's non-GMP BigInteger fallback does not compile at the pinned commit.
find_path(GMP_INCLUDE_DIR gmpxx.h REQUIRED)
find_library(GMP_LIBRARY gmp REQUIRED)
find_library(GMPXX_LIBRARY gmpxx REQUIRED)

# Robustness patch to the pinned Directional (2573873), applied to a build-local copy of the header
# (the submodule stays untouched). In NFunctionMesher::segment_arrangement, coincident arrangement
# edges are meant to be merged whatever their orientation, but the second clause of the test repeats
# the first, so an isoline running exactly along a triangle edge in the opposite direction becomes a
# duplicate halfedge and the triangle's DCEL is inconsistent (empty mesh in Release builds). This
# happens whenever two adjacent singularities are rounded to lattice points that share a coordinate.
set(OVERLAY ${CMAKE_BINARY_DIR}/overlay)
file(READ ${DIRECTIONAL_ROOT}/include/directional/generate_mesh.h GENMESH)
set(BUGGY "((arrEdges[i].first == arrEdges[j].first) && (arrEdges[i].second == arrEdges[j].second))) {")
set(FIXED "((arrEdges[i].first == arrEdges[j].second) && (arrEdges[i].second == arrEdges[j].first))) {")
string(FIND "${GENMESH}" "${BUGGY}" POS)
if(POS EQUAL -1)
  message(FATAL_ERROR "generate_mesh.h changed upstream; re-check the reversed-edge merge patch")
endif()
string(REPLACE "${BUGGY}" "${FIXED}" GENMESH "${GENMESH}")
file(WRITE ${OVERLAY}/directional/generate_mesh.h "${GENMESH}")

add_executable(phcpf_mesher main.cpp)
target_include_directories(phcpf_mesher BEFORE PRIVATE ${OVERLAY})
target_include_directories(phcpf_mesher PRIVATE
  ${DIRECTIONAL_ROOT}/include
  ${DIRECTIONAL_ROOT}/external/eigen
  ${GMP_INCLUDE_DIR})
target_compile_definitions(phcpf_mesher PRIVATE USE_GMP_ENABLED)
target_compile_options(phcpf_mesher PRIVATE -w)
target_link_libraries(phcpf_mesher PRIVATE ${GMPXX_LIBRARY} ${GMP_LIBRARY})
