cmake_minimum_required(VERSION 3.5)

enable_testing()

add_subdirectory(Fem1d)
add_subdirectory(Utils)

add_executable(main main.cpp)
target_link_libraries(main Fem1d)
