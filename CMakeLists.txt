cmake_minimum_required(VERSION 3.15)
project(ToricCodeMWPM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Main executable
add_executable(main main.cpp)

# Test executable
add_executable(test test_main_numerico.cpp)
