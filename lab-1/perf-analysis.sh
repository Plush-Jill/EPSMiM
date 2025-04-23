#!/bin/bash

if [ -z "$1" ]; then
    echo "Ошибка: укажите путь к исполняемому файлу как первый аргумент."
    exit 1
fi

EXECUTABLE="$1"


perf record -g -e cycles -e cache-misses -e cache-references -e branches -e branch-misses $EXECUTABLE
perf report

