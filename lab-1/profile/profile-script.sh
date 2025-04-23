#!/bin/bash

if [ -z "$1" ]; then
    echo "Ошибка: укажите путь к исполняемому файлу как первый аргумент."
    exit 1
fi

EXECUTABLE="$1"

#заменить на билд через CMake
#g++ ../src/main.cpp ../src/heat_source_circle.hpp ../src/heat_source_circle.cpp -funroll-loops -std=c++20 -lboost_json -pg -o ../cmake-build-debug/lab-1
sudo perf record -F 99 -g -- $EXECUTABLE	
sudo perf script > out.perf
../../../FlameGraph/stackcollapse-perf.pl out.perf > out.folded
../../../FlameGraph/flamegraph.pl out.folded  > flamegraph.svg
xdg-open flamegraph.svg




