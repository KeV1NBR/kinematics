# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ical/kevin/library/kinematics

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ical/kevin/library/kinematics/build

# Include any dependencies generated for this target.
include lib/CMakeFiles/kinematics.dir/depend.make

# Include the progress variables for this target.
include lib/CMakeFiles/kinematics.dir/progress.make

# Include the compile flags for this target's objects.
include lib/CMakeFiles/kinematics.dir/flags.make

lib/CMakeFiles/kinematics.dir/kinematics.cpp.o: lib/CMakeFiles/kinematics.dir/flags.make
lib/CMakeFiles/kinematics.dir/kinematics.cpp.o: ../lib/kinematics.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ical/kevin/library/kinematics/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lib/CMakeFiles/kinematics.dir/kinematics.cpp.o"
	cd /home/ical/kevin/library/kinematics/build/lib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/kinematics.dir/kinematics.cpp.o -c /home/ical/kevin/library/kinematics/lib/kinematics.cpp

lib/CMakeFiles/kinematics.dir/kinematics.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/kinematics.dir/kinematics.cpp.i"
	cd /home/ical/kevin/library/kinematics/build/lib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ical/kevin/library/kinematics/lib/kinematics.cpp > CMakeFiles/kinematics.dir/kinematics.cpp.i

lib/CMakeFiles/kinematics.dir/kinematics.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/kinematics.dir/kinematics.cpp.s"
	cd /home/ical/kevin/library/kinematics/build/lib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ical/kevin/library/kinematics/lib/kinematics.cpp -o CMakeFiles/kinematics.dir/kinematics.cpp.s

lib/CMakeFiles/kinematics.dir/kinematics.cpp.o.requires:

.PHONY : lib/CMakeFiles/kinematics.dir/kinematics.cpp.o.requires

lib/CMakeFiles/kinematics.dir/kinematics.cpp.o.provides: lib/CMakeFiles/kinematics.dir/kinematics.cpp.o.requires
	$(MAKE) -f lib/CMakeFiles/kinematics.dir/build.make lib/CMakeFiles/kinematics.dir/kinematics.cpp.o.provides.build
.PHONY : lib/CMakeFiles/kinematics.dir/kinematics.cpp.o.provides

lib/CMakeFiles/kinematics.dir/kinematics.cpp.o.provides.build: lib/CMakeFiles/kinematics.dir/kinematics.cpp.o


# Object files for target kinematics
kinematics_OBJECTS = \
"CMakeFiles/kinematics.dir/kinematics.cpp.o"

# External object files for target kinematics
kinematics_EXTERNAL_OBJECTS =

lib/libkinematics.a: lib/CMakeFiles/kinematics.dir/kinematics.cpp.o
lib/libkinematics.a: lib/CMakeFiles/kinematics.dir/build.make
lib/libkinematics.a: lib/CMakeFiles/kinematics.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ical/kevin/library/kinematics/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libkinematics.a"
	cd /home/ical/kevin/library/kinematics/build/lib && $(CMAKE_COMMAND) -P CMakeFiles/kinematics.dir/cmake_clean_target.cmake
	cd /home/ical/kevin/library/kinematics/build/lib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/kinematics.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/CMakeFiles/kinematics.dir/build: lib/libkinematics.a

.PHONY : lib/CMakeFiles/kinematics.dir/build

lib/CMakeFiles/kinematics.dir/requires: lib/CMakeFiles/kinematics.dir/kinematics.cpp.o.requires

.PHONY : lib/CMakeFiles/kinematics.dir/requires

lib/CMakeFiles/kinematics.dir/clean:
	cd /home/ical/kevin/library/kinematics/build/lib && $(CMAKE_COMMAND) -P CMakeFiles/kinematics.dir/cmake_clean.cmake
.PHONY : lib/CMakeFiles/kinematics.dir/clean

lib/CMakeFiles/kinematics.dir/depend:
	cd /home/ical/kevin/library/kinematics/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ical/kevin/library/kinematics /home/ical/kevin/library/kinematics/lib /home/ical/kevin/library/kinematics/build /home/ical/kevin/library/kinematics/build/lib /home/ical/kevin/library/kinematics/build/lib/CMakeFiles/kinematics.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/CMakeFiles/kinematics.dir/depend
