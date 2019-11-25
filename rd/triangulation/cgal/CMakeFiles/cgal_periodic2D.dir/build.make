# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/benmorton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/triangulation/cgal"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/benmorton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/triangulation/cgal"

# Include any dependencies generated for this target.
include CMakeFiles/cgal_periodic2D.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/cgal_periodic2D.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cgal_periodic2D.dir/flags.make

CMakeFiles/cgal_periodic2D.dir/cgal_periodic2D.cpp.o: CMakeFiles/cgal_periodic2D.dir/flags.make
CMakeFiles/cgal_periodic2D.dir/cgal_periodic2D.cpp.o: cgal_periodic2D.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/benmorton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/triangulation/cgal/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/cgal_periodic2D.dir/cgal_periodic2D.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cgal_periodic2D.dir/cgal_periodic2D.cpp.o -c "/Users/benmorton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/triangulation/cgal/cgal_periodic2D.cpp"

CMakeFiles/cgal_periodic2D.dir/cgal_periodic2D.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cgal_periodic2D.dir/cgal_periodic2D.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/benmorton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/triangulation/cgal/cgal_periodic2D.cpp" > CMakeFiles/cgal_periodic2D.dir/cgal_periodic2D.cpp.i

CMakeFiles/cgal_periodic2D.dir/cgal_periodic2D.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cgal_periodic2D.dir/cgal_periodic2D.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/benmorton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/triangulation/cgal/cgal_periodic2D.cpp" -o CMakeFiles/cgal_periodic2D.dir/cgal_periodic2D.cpp.s

# Object files for target cgal_periodic2D
cgal_periodic2D_OBJECTS = \
"CMakeFiles/cgal_periodic2D.dir/cgal_periodic2D.cpp.o"

# External object files for target cgal_periodic2D
cgal_periodic2D_EXTERNAL_OBJECTS =

cgal_periodic2D: CMakeFiles/cgal_periodic2D.dir/cgal_periodic2D.cpp.o
cgal_periodic2D: CMakeFiles/cgal_periodic2D.dir/build.make
cgal_periodic2D: /opt/local/lib/libCGAL.13.0.2.dylib
cgal_periodic2D: /opt/local/lib/libgmp.dylib
cgal_periodic2D: /opt/local/lib/libmpfr.dylib
cgal_periodic2D: /opt/local/lib/libboost_thread-mt.dylib
cgal_periodic2D: /opt/local/lib/libboost_system-mt.dylib
cgal_periodic2D: /opt/local/lib/libboost_chrono-mt.dylib
cgal_periodic2D: /opt/local/lib/libboost_date_time-mt.dylib
cgal_periodic2D: /opt/local/lib/libboost_atomic-mt.dylib
cgal_periodic2D: /usr/local/lib/libboost_thread-mt.dylib
cgal_periodic2D: CMakeFiles/cgal_periodic2D.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/benmorton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/triangulation/cgal/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable cgal_periodic2D"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cgal_periodic2D.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cgal_periodic2D.dir/build: cgal_periodic2D

.PHONY : CMakeFiles/cgal_periodic2D.dir/build

CMakeFiles/cgal_periodic2D.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cgal_periodic2D.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cgal_periodic2D.dir/clean

CMakeFiles/cgal_periodic2D.dir/depend:
	cd "/Users/benmorton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/triangulation/cgal" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/benmorton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/triangulation/cgal" "/Users/benmorton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/triangulation/cgal" "/Users/benmorton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/triangulation/cgal" "/Users/benmorton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/triangulation/cgal" "/Users/benmorton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/triangulation/cgal/CMakeFiles/cgal_periodic2D.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/cgal_periodic2D.dir/depend

