#!/bin/bash
#
# ./cmake-build-debug/lab-1 это я исполнял скрипт из директории на 1 выше директории с бинарником, заменишь на своё.
/opt/intel/oneapi/advisor/2025.0/bin64/advixe-cl --collect=survey --project-dir=./ -- ./cmake-build-debug/lab-1
/opt/intel/oneapi/advisor/2025.0/bin64/advixe-cl --collect=tripcounts --flop --project-dir=./ -- ./cmake-build-debug/lab-1
/opt/intel/oneapi/advisor/2025.0/bin64/advixe-cl --report=roofline --project-dir=./
