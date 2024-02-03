# test/CMakeLists.txt

# Add test source files in the test directory
file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# Link with the main library
target_link_libraries(mddriver_lib)

# Create the test executable
add_executable(test_mddriver ${TEST_SOURCES})

# Link with the test executable
target_link_libraries(test_mddriver mddriver_lib)

