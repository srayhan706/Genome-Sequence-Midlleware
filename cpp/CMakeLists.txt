cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)

find_package(MPI REQUIRED COMPONENTS CXX)
find_package(OpenMP)
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

add_executable(genomic_middleware
  src/main.cpp
  src/server.cpp
  src/worker.cpp
  src/dispatch.cpp
  src/dispatch.hpp
  src/mpi_protocol.cpp
  src/mpi_protocol.hpp
  src/algorithms/smith_waterman.cpp
  src/algorithms/smith_waterman.hpp
  src/algorithms/kmp.cpp
  src/algorithms/kmp.hpp
  src/algorithms/kmer.cpp
  src/algorithms/kmer.hpp
)

target_include_directories(genomic_middleware PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${GEN_PATH}
  ${Protobuf_INCLUDE_DIRS}
)
target_link_libraries(genomic_middleware PRIVATE genomic_proto MPI::MPI_CXX)
if(OpenMP_CXX_FOUND)
  target_link_libraries(genomic_middleware PRIVATE OpenMP::OpenMP_CXX)
  target_compile_definitions(genomic_middleware PRIVATE _OPENMP)
endif()
