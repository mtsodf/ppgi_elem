add_executable(LinearAlgebraTest LinearAlgebraTest.cpp)

target_link_libraries(LinearAlgebraTest Solvers Utils)

add_test(NAME LinearAlgebraTest COMMAND "${CMAKE_CURRENT_BINARY_DIR}/LinearAlgebraTest")