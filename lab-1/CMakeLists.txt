cmake_minimum_required(VERSION 3.30)
project(lab-1)

set(CMAKE_CXX_STANDARD 20)


add_executable(lab-1
        src/common/main.cpp
        src/common/heat_source_circle.cpp
        src/common/heat_source_circle.hpp
        src/lab-1/poisson_equation_solver.cpp
        src/common/poisson_equation_solver.hpp
)
target_compile_options(lab-1 PRIVATE -O3)
target_link_libraries(lab-1
    boost_json
)
target_compile_definitions(lab-1 PRIVATE LAB_1)


add_executable(lab-2
        src/common/main.cpp
        src/common/heat_source_circle.cpp
        src/common/heat_source_circle.hpp
        src/lab-2/poisson_equation_solver.cpp
        src/lab-2/poisson_equation_solver.hpp
        src/lab-2/aligned_allocator.hpp
)

#target_compile_options(lab-2 PRIVATE -O3 -fopenmp -ftree-vectorize -march=native -fopt-info-vec -fopenmp-simd -mavx512f)
target_compile_options(lab-2 PRIVATE -O3 -fopenmp -mavx512f -ftree-vectorize -march=native -fopenmp-simd)
#target_compile_options(lab-2 PRIVATE -DNSIMD_AVX2 -O3 -I/home/plush-jill/All_Random/git/nsimd/include -march=native)
find_package(OpenMP REQUIRED)
target_link_libraries(lab-2
        boost_json
        OpenMP::OpenMP_CXX)
target_compile_definitions(lab-2 PRIVATE LAB_2)

#if(OpenMP_CXX_FOUND)
#    target_link_libraries(lab-2 PUBLIC OpenMP::OpenMP_CXX)
#endif()


add_executable(stuff
        stuff/main.cpp
)
