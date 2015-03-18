cmake_minimum_required(VERSION 2.8)
project(OpenTspline)

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include/${PROJECT_NAME} CACHE PATH
  "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/cmake/${PROJECT_NAME})
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

add_definitions("-Wall -Wno-unused-function -fno-strict-aliasing -pthread")

include_directories(.)

include(openmp.cmake)

set(LIBRARY_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib")

add_subdirectory(Tspline)
add_subdirectory(apps)


# OpenTsplineConfig.cmake
# -----------------------
# Add all targets to the build-tree export set
export(TARGETS tspline FILE "${PROJECT_BINARY_DIR}/OpenTsplineTargets.cmake")
 
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE OpenTspline)
 
# Create the FooBarConfig.cmake and FooBarConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}")
set(CONF_LIBRARY_DIRS "${PROJECT_SOURCE_DIR}/lib")
configure_file(OpenTsplineConfig.cmake.in "${PROJECT_BINARY_DIR}/OpenTsplineConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "${INSTALL_INCLUDE_DIR}")
set(CONF_LIBRARY_DIRS "${INSTALL_LIB_DIR}")
configure_file(OpenTsplineConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/OpenTsplineConfig.cmake" @ONLY)
# ... for both
configure_file(OpenTsplineConfigVersion.cmake.in "${PROJECT_BINARY_DIR}/OpenTsplineConfigVersion.cmake" @ONLY)
 
# Install the OpenTsplineConfig.cmake and OpenTsplineConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/OpenTsplineConfig.cmake"
  "${PROJECT_BINARY_DIR}/OpenTsplineConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
 
# Install the export set for use with the install-tree
install(EXPORT OpenTsplineTargets DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)



# uninstall target
# ----------------
configure_file(
"${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
"${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
IMMEDIATE @ONLY)
add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
