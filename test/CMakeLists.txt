# test/CMakeLists.txt

# Add test source files in the test directory
# file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
find_package(Eigen3 REQUIRED)

# Link with the main library
target_link_libraries(mddriver_lib)

# Create the test executable
add_executable(rk45-lorenz "rk45-lorenz.cpp")

# Link with the test executable
target_link_libraries(rk45-lorenz PRIVATE
    mddriver_lib 
    Eigen3::Eigen
)

# Create the test executable
add_executable(rk4-lorenz "rk4-lorenz.cpp")

# Link with the test executable
target_link_libraries(rk4-lorenz PRIVATE
    mddriver_lib 
    Eigen3::Eigen
)

