#!/bin/bash


g++ ../src/main.cpp ../src/heat_source_circle.hpp ../src/heat_source_circle.cpp -funroll-loops -std=c++20 -lboost_json -pg -o ../cmake-build-debug/lab-1
sudo perf record -F 99 -g -- ../cmake-build-debug/lab-1	
sudo perf script > out.perf
../../../FlameGraph/stackcollapse-perf.pl out.perf > out.folded
../../../FlameGraph/flamegraph.pl out.folded  > flamegraph.svg
xdg-open flamegraph.svg



$$
\begin{pmatrix}
1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 1 & 0 \\
0 & 0 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 1 & 0 & 0 & 1 & 0 & 0 & 1 \\
0 & 0 & 1 & 0 & 1 & 1 & 0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 1 & 0 & 1 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 0 & 0 & 0 \\
\end{pmatrix}
$$