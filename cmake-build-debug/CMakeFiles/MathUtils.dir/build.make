# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

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
CMAKE_COMMAND = /root/clion-2020.2.1/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /root/clion-2020.2.1/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /root/Desktop/MathUtils

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /root/Desktop/MathUtils/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/MathUtils.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/MathUtils.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MathUtils.dir/flags.make

CMakeFiles/MathUtils.dir/main.cpp.o: CMakeFiles/MathUtils.dir/flags.make
CMakeFiles/MathUtils.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/root/Desktop/MathUtils/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MathUtils.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MathUtils.dir/main.cpp.o -c /root/Desktop/MathUtils/main.cpp

CMakeFiles/MathUtils.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MathUtils.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /root/Desktop/MathUtils/main.cpp > CMakeFiles/MathUtils.dir/main.cpp.i

CMakeFiles/MathUtils.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MathUtils.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /root/Desktop/MathUtils/main.cpp -o CMakeFiles/MathUtils.dir/main.cpp.s

CMakeFiles/MathUtils.dir/library.cpp.o: CMakeFiles/MathUtils.dir/flags.make
CMakeFiles/MathUtils.dir/library.cpp.o: ../library.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/root/Desktop/MathUtils/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/MathUtils.dir/library.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MathUtils.dir/library.cpp.o -c /root/Desktop/MathUtils/library.cpp

CMakeFiles/MathUtils.dir/library.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MathUtils.dir/library.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /root/Desktop/MathUtils/library.cpp > CMakeFiles/MathUtils.dir/library.cpp.i

CMakeFiles/MathUtils.dir/library.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MathUtils.dir/library.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /root/Desktop/MathUtils/library.cpp -o CMakeFiles/MathUtils.dir/library.cpp.s

CMakeFiles/MathUtils.dir/DataType.cpp.o: CMakeFiles/MathUtils.dir/flags.make
CMakeFiles/MathUtils.dir/DataType.cpp.o: ../DataType.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/root/Desktop/MathUtils/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/MathUtils.dir/DataType.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MathUtils.dir/DataType.cpp.o -c /root/Desktop/MathUtils/DataType.cpp

CMakeFiles/MathUtils.dir/DataType.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MathUtils.dir/DataType.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /root/Desktop/MathUtils/DataType.cpp > CMakeFiles/MathUtils.dir/DataType.cpp.i

CMakeFiles/MathUtils.dir/DataType.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MathUtils.dir/DataType.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /root/Desktop/MathUtils/DataType.cpp -o CMakeFiles/MathUtils.dir/DataType.cpp.s

# Object files for target MathUtils
MathUtils_OBJECTS = \
"CMakeFiles/MathUtils.dir/main.cpp.o" \
"CMakeFiles/MathUtils.dir/library.cpp.o" \
"CMakeFiles/MathUtils.dir/DataType.cpp.o"

# External object files for target MathUtils
MathUtils_EXTERNAL_OBJECTS =

MathUtils: CMakeFiles/MathUtils.dir/main.cpp.o
MathUtils: CMakeFiles/MathUtils.dir/library.cpp.o
MathUtils: CMakeFiles/MathUtils.dir/DataType.cpp.o
MathUtils: CMakeFiles/MathUtils.dir/build.make
MathUtils: CMakeFiles/MathUtils.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/root/Desktop/MathUtils/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable MathUtils"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MathUtils.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MathUtils.dir/build: MathUtils

.PHONY : CMakeFiles/MathUtils.dir/build

CMakeFiles/MathUtils.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MathUtils.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MathUtils.dir/clean

CMakeFiles/MathUtils.dir/depend:
	cd /root/Desktop/MathUtils/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /root/Desktop/MathUtils /root/Desktop/MathUtils /root/Desktop/MathUtils/cmake-build-debug /root/Desktop/MathUtils/cmake-build-debug /root/Desktop/MathUtils/cmake-build-debug/CMakeFiles/MathUtils.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/MathUtils.dir/depend

