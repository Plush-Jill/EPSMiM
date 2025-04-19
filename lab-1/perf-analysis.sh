#!/bin/bash

perf record -g -e cycles -e cache-misses -e cache-references -e branches -e branch-misses ./cmake-build-debug/lab-1
perf report

