cmake_minimum_required(VERSION 2.8)

SET(LIBRARY_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib")
SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/bin")

add_definitions("-Wall -Wno-unused-function -fno-strict-aliasing -pthread")

include_directories(.)

include(openmp.cmake)

add_subdirectory(v4r)
add_subdirectory(apps)

# uninstall target
configure_file(
"${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
"${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
IMMEDIATE @ONLY)
add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
