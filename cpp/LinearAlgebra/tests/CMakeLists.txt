add_executable(LinearAlgebraTest LinearAlgebraTest.cpp)

target_link_libraries(LinearAlgebraTest Solvers Utils)

add_test(NAME LinearAlgebraTest COMMAND "${CMAKE_CURRENT_BINARY_DIR}/LinearAlgebraTest")



include_directories(${CMAKE_SOURCE_DIR}/ThirdParty/viennacl)
include_directories(../../ThirdParty/)
add_executable(ViennaClTest ViennaClTest.cpp)
#add_executable
#find_package(CUDA QUIET REQUIRED)
add_executable(sparse sparse.cpp)
configure_file(mat65k.mtx mat65k.mtx COPYONLY)