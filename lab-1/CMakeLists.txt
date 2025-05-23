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
        src/common/aligned_allocator.hpp
)

#target_compile_options(lab-2 PRIVATE -O3 -fopenmp -ftree-vectorize -march=native -fopt-info-vec -fopenmp-simd -mavx512f)
target_compile_options(lab-2 PRIVATE -O3 -fopenmp -mavx512f -ftree-vectorize -march=native -fopenmp-simd -funroll-loops)
#target_compile_options(lab-2 PRIVATE -DNSIMD_AVX2 -O3 -I/home/plush-jill/All_Random/git/nsimd/include -march=native)
find_package(OpenMP REQUIRED)
target_link_libraries(lab-2
        boost_json
        OpenMP::OpenMP_CXX)
target_compile_definitions(lab-2 PRIVATE LAB_2)


add_executable(lab-3
        src/common/main.cpp
        src/common/heat_source_circle.cpp
        src/common/heat_source_circle.hpp
        src/common/aligned_allocator.hpp
        src/lab-3/poisson_equation_solver.cpp
        src/lab-3/poisson_equation_solver.hpp
        src/lab-3/front_moving_along_array.cpp
        src/lab-3/front_moving_along_array.hpp
)

#target_compile_options(lab-2 PRIVATE -O3 -fopenmp -ftree-vectorize -march=native -fopt-info-vec -fopenmp-simd -mavx512f)
target_compile_options(lab-3 PRIVATE -O3 -fopenmp -mavx512f -ftree-vectorize -march=native -fopenmp-simd)
#target_compile_options(lab-3 PRIVATE -DNSIMD_AVX2 -O3 -I/home/plush-jill/All_Random/git/nsimd/include -march=native)
find_package(OpenMP REQUIRED)
target_link_libraries(lab-3
        boost_json
        OpenMP::OpenMP_CXX)
target_compile_definitions(lab-3 PRIVATE LAB_3)



add_executable(lab-4
        src/common/main.cpp
        src/common/heat_source_circle.cpp
        src/common/heat_source_circle.hpp
        src/common/aligned_allocator.hpp
        src/lab-4/poisson_equation_solver.cpp
        src/lab-4/poisson_equation_solver.hpp
        src/lab-4/window/left/window_left_part.cpp
        src/lab-4/window/left/window_left_part.hpp
        src/lab-4/window/window_abstract.cpp
        src/lab-4/window/window_abstract.hpp
        src/lab-4/window/middle/window_middle_part.cpp
        src/lab-4/window/middle/window_middle_part.hpp
        src/lab-4/window/right/window_right_part.cpp
        src/lab-4/window/right/window_right_part.hpp
)

target_compile_options(lab-4 PRIVATE -O3 -fopenmp -mavx512f -ftree-vectorize -march=native -fopenmp-simd)
find_package(OpenMP REQUIRED)
target_link_libraries(lab-4
        boost_json
        OpenMP::OpenMP_CXX)
target_compile_definitions(lab-4 PRIVATE LAB_4)


#if(OpenMP_CXX_FOUND)
#    target_link_libraries(lab-2 PUBLIC OpenMP::OpenMP_CXX)
#endif()


add_executable(stuff
        stuff/stuff.cpp
        src/common/aligned_allocator.hpp
)
