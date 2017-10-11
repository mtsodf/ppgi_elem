cmake_minimum_required(VERSION 3.5)

enable_testing()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(OpenCL REQUIRED)

message("OpenCL include dir: ${OPENCL_INCLUDE_DIR}")
message("OpenCL library dir: ${OPENCL_LIBRARY}")

add_subdirectory(Fem1d)
add_subdirectory(LinearAlgebra)
add_subdirectory(Utils)

add_executable(main main.cpp)
target_link_libraries(main Fem1d)
