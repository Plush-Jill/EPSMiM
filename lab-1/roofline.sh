#!/bin/bash

if [ -z "$1" ]; then
    echo "Ошибка: укажите путь к исполняемому файлу как первый аргумент."
    exit 1
fi

EXECUTABLE="$1"
/opt/intel/oneapi/advisor/2025.1/bin64/advixe-cl --collect=survey --project-dir=./ -- $EXECUTABLE
/opt/intel/oneapi/advisor/2025.1/bin64/advixe-cl --collect=tripcounts --flop --project-dir=./ -- $EXECUTABLE
/opt/intel/oneapi/advisor/2025.1/bin64/advixe-cl --report=roofline --project-dir=./

