configure_file(PlotFem1d.py PlotFem1d.py COPYONLY)

add_library(Fem1d SHARED Fem1d.cpp )

target_link_libraries(Fem1d Solvers Utils)
add_subdirectory(tests)