cmake_minimum_required(VERSION 3.10)

project(sem10_part1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(app
    src/main.cpp
    src/methods.cpp
    src/csv_output.cpp
    src/benchmark.cpp          # <-- добавили
)

target_include_directories(app PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
