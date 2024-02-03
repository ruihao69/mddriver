cmake_minimum_required(VERSION 3.12)

project(mddriver)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Find dependencies
# find_package(Eigen3 REQUIRED)

# Set include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
)

# Add source files in the src directory
# file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")

# Create the executable
# add_executable(mddriver ${SOURCES})

# Add subdirectories for specific components
add_subdirectory(src)
add_subdirectory(test)

