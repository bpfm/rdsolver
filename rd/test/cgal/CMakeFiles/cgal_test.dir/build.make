# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.13.4/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.13.4/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/morton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/test/cgal"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/morton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/test/cgal"

# Include any dependencies generated for this target.
include CMakeFiles/cgal_test.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/cgal_test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cgal_test.dir/flags.make

CMakeFiles/cgal_test.dir/cgal_test.cpp.o: CMakeFiles/cgal_test.dir/flags.make
CMakeFiles/cgal_test.dir/cgal_test.cpp.o: cgal_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/morton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/test/cgal/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/cgal_test.dir/cgal_test.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cgal_test.dir/cgal_test.cpp.o -c "/Users/morton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/test/cgal/cgal_test.cpp"

CMakeFiles/cgal_test.dir/cgal_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cgal_test.dir/cgal_test.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/morton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/test/cgal/cgal_test.cpp" > CMakeFiles/cgal_test.dir/cgal_test.cpp.i

CMakeFiles/cgal_test.dir/cgal_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cgal_test.dir/cgal_test.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/morton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/test/cgal/cgal_test.cpp" -o CMakeFiles/cgal_test.dir/cgal_test.cpp.s

# Object files for target cgal_test
cgal_test_OBJECTS = \
"CMakeFiles/cgal_test.dir/cgal_test.cpp.o"

# External object files for target cgal_test
cgal_test_EXTERNAL_OBJECTS =

cgal_test: CMakeFiles/cgal_test.dir/cgal_test.cpp.o
cgal_test: CMakeFiles/cgal_test.dir/build.make
cgal_test: /usr/local/lib/libCGAL.13.0.3.dylib
cgal_test: /usr/local/lib/libmpfr.dylib
cgal_test: /usr/local/lib/libgmp.dylib
cgal_test: /usr/local/lib/libboost_thread-mt.dylib
cgal_test: /usr/local/lib/libboost_system-mt.dylib
cgal_test: /usr/local/lib/libboost_chrono-mt.dylib
cgal_test: /usr/local/lib/libboost_date_time-mt.dylib
cgal_test: /usr/local/lib/libboost_atomic-mt.dylib
cgal_test: /usr/local/lib/libboost_thread-mt.dylib
cgal_test: /usr/local/lib/libboost_chrono-mt.dylib
cgal_test: /usr/local/lib/libboost_date_time-mt.dylib
cgal_test: /usr/local/lib/libboost_atomic-mt.dylib
cgal_test: CMakeFiles/cgal_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/morton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/test/cgal/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable cgal_test"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cgal_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cgal_test.dir/build: cgal_test

.PHONY : CMakeFiles/cgal_test.dir/build

CMakeFiles/cgal_test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cgal_test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cgal_test.dir/clean

CMakeFiles/cgal_test.dir/depend:
	cd "/Users/morton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/test/cgal" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/morton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/test/cgal" "/Users/morton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/test/cgal" "/Users/morton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/test/cgal" "/Users/morton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/test/cgal" "/Users/morton/Library/Mobile Documents/com~apple~CloudDocs/rdsolver/rd/test/cgal/CMakeFiles/cgal_test.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/cgal_test.dir/depend

