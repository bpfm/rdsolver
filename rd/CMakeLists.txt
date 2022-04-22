cmake_minimum_required(VERSION 3.17)
project(MovingMeshRD)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
set(CMAKE_BUILD_TYPE Release)

include_directories(.)
include_directories("/usr/local/Cellar/openblas/0.3.19/include")
include_directories("/usr/local/Cellar/libomp/13.0.0/include")
set(OpenBLAS_DIR "/usr/local/Cellar/openblas/0.3.19/lib/cmake/openblas")


find_package(OpenBLAS REQUIRED)
find_package(OpenMP REQUIRED)

#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()



# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )
if ( NOT CGAL_FOUND )
    message(STATUS "This project requires the CGAL library, and will not be compiled.")
    return()
endif()


# Boost and its components
find_package( Boost REQUIRED )
if ( NOT Boost_FOUND )
    message(STATUS "This project requires the Boost library, and will not be compiled.")
    return()
endif()

add_executable(MovingMeshRD
        src/base.cpp
        src/constants.h
        src/inverse.cpp
        src/io2D.cpp
        src/main.cpp
        src/setup2D.cpp
        src/source2D.cpp
        src/timestep.cpp
        src/triangle2D.h
        src/vertex2D.h
        src/all_functions.h
        src/base.h
        src/inverse.h)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS MovingMeshRD )

target_link_libraries(MovingMeshRD PRIVATE CGAL::CGAL ${OpenBLAS_LIBRARIES} OpenMP::OpenMP_CXX)

