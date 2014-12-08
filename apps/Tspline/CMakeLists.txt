set(APPS_TSPLINE 0 CACHE BOOL "Examples for Tspline library (requires TomGine)" )
IF(APPS_TSPLINE)

# CGAL "Computational Geometry Algorithms Library" dependent programs
find_package(CGAL QUIET COMPONENTS Core )
IF(CGAL_FOUND)

    add_definitions(-DGL_GLEXT_PROTOTYPES)

    # find Eigen
    find_package(PkgConfig REQUIRED)
    pkg_search_module(EIGEN3 REQUIRED eigen3)
    include_directories(${EIGEN3_INCLUDE_DIRS})

    # find OpenCV
    find_package( OpenCV REQUIRED )

    add_executable(TsplineCox Cox.cpp)
    target_link_libraries(TsplineCox tspline v4rTomGine)

    add_executable(Tspline Tspline.cpp)
    target_link_libraries(Tspline tspline v4rTomGine ${OpenCV_LIBRARIES})

    add_executable(TsplineBox TsplineBox.cpp)
    target_link_libraries(TsplineBox tspline v4rTomGine ${OpenCV_LIBRARIES})

    add_executable(TsplineCylinder TsplineCylinder.cpp)
    target_link_libraries(TsplineCylinder tspline v4rTomGine ${OpenCV_LIBRARIES})

    add_executable(TsplineMP2Tsplines TsplineMP2Tsplines.cpp)
    target_link_libraries(TsplineMP2Tsplines tspline v4rTomGine ${OpenCV_LIBRARIES})

    add_executable(TsplineRefinement TsplineRefinement.cpp)
    target_link_libraries(TsplineRefinement tspline v4rTomGine ${OpenCV_LIBRARIES})

    add_executable(TsplineSequenceViewer TsplineSequenceViewer.cpp)
    target_link_libraries(TsplineSequenceViewer tspline v4rTomGine ${OpenCV_LIBRARIES})

    add_executable(TsplineViewer TsplineViewer.cpp)
    target_link_libraries(TsplineViewer tspline v4rTomGine ${OpenCV_LIBRARIES})

    add_executable(TsplineViewerMP TsplineViewerMP.cpp)
    target_link_libraries(TsplineViewerMP tspline v4rTomGine ${OpenCV_LIBRARIES})

ELSE(CGAL_FOUND)
    message(STATUS "APPS_TSPLINE: requires the CGAL library, and will not be compiled.")
ENDIF(CGAL_FOUND)

ENDIF(APPS_TSPLINE)
