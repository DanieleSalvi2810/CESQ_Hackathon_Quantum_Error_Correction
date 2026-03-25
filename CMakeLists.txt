cmake_minimum_required(VERSION 3.15)
project(ToricCodeMWPM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Main executable
add_executable(main main.cpp)

<<<<<<< Updated upstream
=======

>>>>>>> Stashed changes
# Convenience target: build and run the main executable.
add_custom_target(run-main
    COMMAND $<TARGET_FILE:main>
    DEPENDS main
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running main executable"
)
