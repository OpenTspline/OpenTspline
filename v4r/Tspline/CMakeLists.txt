set(V4R_TSPLINE 1 CACHE BOOL "T-spline implementation" )
IF(V4R_TSPLINE)
project(v4rTspline)

# CGAL "Computational Geometry Algorithms Library" dependent programs
find_package(CGAL QUIET COMPONENTS Core )

set(CMAKE_CXX_FLAGS "-frounding-math")

if ( CGAL_FOUND )

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
	set(CMAKE_BUILD_TYPE Release)
	
	# find dependent packages
	find_package(PkgConfig REQUIRED)
	pkg_search_module(EIGEN3 REQUIRED eigen3)
	include_directories(${EIGEN3_INCLUDE_DIRS})

  # create tspline include directory name
  string(REGEX REPLACE "v4r" "" TSPLINE_INCLUDE_NAME ${PROJECT_NAME})

	include( ${CGAL_USE_FILE} )
	
	set(SOURCE_CPP
                Cox.cpp
                Domain.cpp
                File.cpp
                Polyhedron.cpp
                Transition.cpp
                Tspline.cpp
                TsplineCreator.cpp
                Tspline_operators.cpp
                TsplineMultiPatch.cpp
                TsplinePatch.cpp
                Types.cpp)

	set(SOURCE_H
                Cox.h
                Domain.h
                File.h
                Given_border_parameterizer_3.h
                Math.hpp
                Polyhedron.h
                Transition.h
                Tspline.h
                TsplineCreator.h
                TsplineMultiPatch.h
                TsplinePatch.h
                Types.h)

  add_library(${PROJECT_NAME} SHARED ${SOURCE_CPP} ${SOURCE_H})
  target_link_libraries(${PROJECT_NAME} ${CGAL_LIBRARIES})


  if(${V4R_BUILD})
    v4r_add_library(${PROJECT_NAME} "${SOURCE_H}")
  else(${V4R_BUILD})
    install(DIRECTORY DESTINATION include/v4r/${TSPLINE_INCLUDE_NAME})
    install(FILES ${SOURCE_H} DESTINATION include/v4r/${TSPLINE_INCLUDE_NAME})
    install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION lib)
  endif(${V4R_BUILD})


else(CGAL_FOUND)
    message(STATUS "V4R_TSPLINE: class Tspline requires the CGAL library, and will not be compiled.")
endif(CGAL_FOUND)


ENDIF(V4R_TSPLINE)
