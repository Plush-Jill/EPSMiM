cmake_minimum_required(VERSION 3.30)
project(lab_1)

set(CMAKE_CXX_STANDARD 20)

add_executable(lab_1 src/main.cpp
        src/heat_source_circle.cpp
        src/heat_source_circle.hpp)

target_link_libraries(lab_1
boost_json)
