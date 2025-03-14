cmake_minimum_required(VERSION 3.30)
project(lab-1)

set(CMAKE_CXX_STANDARD 20)

#find_package(NSIMD)
#find_package(NSIMD COMPONENTS avx512)

add_executable(lab-1 src/main.cpp
        src/heat_source_circle.cpp
        src/heat_source_circle.hpp
        src/poisson_equation_solver.cpp
        src/poisson_equation_solver.hpp)
add_compile_options(-O3)
target_link_libraries(lab-1
boost_json)
