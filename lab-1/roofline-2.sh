#!/bin/bash
/opt/intel/oneapi/advisor/2025.0/bin64/advixe-cl --collect=survey --project-dir=./ -- ./cmake-build-debug/lab-2
/opt/intel/oneapi/advisor/2025.0/bin64/advixe-cl --collect=tripcounts --flop --project-dir=./ -- ./cmake-build-debug/lab-2
/opt/intel/oneapi/advisor/2025.0/bin64/advixe-cl --report=roofline --project-dir=./
