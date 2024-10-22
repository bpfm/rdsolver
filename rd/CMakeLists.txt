cmake_minimum_required(VERSION 3.17)
project(rd)

# Set the compilers
set(CMAKE_C_COMPILER "/usr/local/opt/llvm/bin/clang")
set(CMAKE_CXX_COMPILER "/usr/local/opt/llvm/bin/clang++")

# Set the standard and flags
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

# Include directories
include_directories(.)
include_directories("/usr/include")
include_directories("/usr/local/include")
include_directories("/usr/local/Cellar/openblas/0.3.28/include")
include_directories("/usr/local/Cellar/libomp/13.0.1/include")

# Set the SDK root
set(CMAKE_OSX_SYSROOT "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk")

# Find packages
set(OpenBLAS_DIR "/usr/local/Cellar/openblas/0.3.28/lib/cmake/openblas")
find_package(OpenBLAS REQUIRED)
find_package(OpenMP REQUIRED)

# Link OpenMP
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Add executables and link libraries
# add_executable(rd2D
#     base.cpp
#     constants.h
#     inverse.cpp
#     io2D.cpp
#     main.cpp
#     setup2D.cpp
#     source2D.cpp
#     timestep.cpp
#     triangle2D.h
#     vertex2D.h
#     all_functions.h
#     base.h
#     inverse.h
# )
# target_link_libraries(rd2D ${OpenBLAS_LIBRARIES} OpenMP::OpenMP_CXX)

add_executable(rd3D
    base.cpp
    constants.h
    inverse.cpp
    io3D.cpp
    main.cpp
    setup3D.cpp
    source3D.cpp
    timestep.cpp
    triangle3D.h
    vertex3D.h
    all_functions.h
    base.h
    inverse.h
)
target_link_libraries(rd3D ${OpenBLAS_LIBRARIES} OpenMP::OpenMP_CXX)
