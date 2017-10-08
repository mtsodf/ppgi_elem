add_executable(Fem1dTest Fem1dTest.cpp)

target_link_libraries(Fem1dTest Fem1d Utils)

add_test(NAME Fem1dTest COMMAND "${CMAKE_CURRENT_BINARY_DIR}/Fem1dTest")