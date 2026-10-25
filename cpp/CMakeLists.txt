cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)

find_package(MPI REQUIRED COMPONENTS CXX)

# Apple Clang has no bundled OpenMP runtime; point FindOpenMP at Homebrew's libomp.
if(APPLE AND NOT OpenMP_ROOT)
  execute_process(COMMAND brew --prefix libomp
    OUTPUT_VARIABLE _libomp_prefix OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(_libomp_prefix AND EXISTS "${_libomp_prefix}")
    set(OpenMP_ROOT "${_libomp_prefix}")
  endif()
endif()
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP enabled (${OpenMP_CXX_FLAGS})")
else()
  message(WARNING "OpenMP not found: algorithms will run single-threaded. On macOS: brew install libomp")
endif()
find_package(Protobuf REQUIRED)
find_program(GRPC_CPP_PLUGIN grpc_cpp_plugin
  HINTS $ENV{HOME}/.local/bin /usr/local/bin /opt/homebrew/bin)

if(NOT GRPC_CPP_PLUGIN)
  message(FATAL_ERROR "grpc_cpp_plugin not found. Install gRPC (e.g. brew install grpc).")
endif()

# Proto generation
set(PROTO_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../proto")
set(GEN_PATH "${CMAKE_CURRENT_BINARY_DIR}/gen")
file(MAKE_DIRECTORY "${GEN_PATH}")

set(GEN_PROTO
  "${GEN_PATH}/genomic.pb.cc" "${GEN_PATH}/genomic.pb.h"
  "${GEN_PATH}/genomic.grpc.pb.cc" "${GEN_PATH}/genomic.grpc.pb.h")

add_custom_command(
  OUTPUT ${GEN_PROTO}
  COMMAND ${Protobuf_PROTOC_EXECUTABLE}
    --cpp_out=${GEN_PATH}
    --grpc_out=${GEN_PATH}
    --plugin=protoc-gen-grpc=${GRPC_CPP_PLUGIN}
    -I${PROTO_PATH}
    ${PROTO_PATH}/genomic.proto
  DEPENDS ${PROTO_PATH}/genomic.proto
  COMMENT "Generating gRPC C++ from genomic.proto"
)

add_library(genomic_proto ${GEN_PROTO})
target_include_directories(genomic_proto PUBLIC ${GEN_PATH} ${Protobuf_INCLUDE_DIRS})
target_link_libraries(genomic_proto PUBLIC protobuf::libprotobuf)

# gRPC libraries via pkg-config (pulls in grpc++, grpc, gpr, abseil)
find_package(PkgConfig REQUIRED)
pkg_check_modules(GRPC REQUIRED grpc++)
target_include_directories(genomic_proto PUBLIC ${GRPC_INCLUDE_DIRS})
target_link_directories(genomic_proto PUBLIC ${GRPC_LIBRARY_DIRS})
target_link_libraries(genomic_proto PUBLIC ${GRPC_LIBRARIES})
target_compile_options(genomic_proto PUBLIC ${GRPC_CFLAGS_OTHER})

# Algorithms + MPI wire format: no MPI/gRPC dependency, so tests can link it alone.
add_library(genomic_core STATIC
  src/mpi_protocol.cpp
  src/algorithms/smith_waterman.cpp
  src/algorithms/kmp.cpp
  src/algorithms/kmer.cpp
)
target_include_directories(genomic_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
if(OpenMP_CXX_FOUND)
  target_link_libraries(genomic_core PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(genomic_middleware
  src/main.cpp
  src/server.cpp
  src/worker.cpp
  src/dispatch.cpp
)
target_link_libraries(genomic_middleware PRIVATE genomic_core genomic_proto MPI::MPI_CXX)

add_executable(algorithms_test tests/algorithms_test.cpp)
target_link_libraries(algorithms_test PRIVATE genomic_core)
add_test(NAME algorithms_test COMMAND algorithms_test)
