# src/CMakeLists.txt

# Add source files in the src directory
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# Add header files in the include directory
include_directories(${CMAKE_SOURCE_DIR}/include)

# Create a library from the source files
add_library(mddriver_lib ${SOURCES})

# Link the library with the main executable
target_link_libraries(mddriver mddriver_lib)

