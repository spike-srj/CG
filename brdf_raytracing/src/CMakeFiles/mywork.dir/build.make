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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/leshui/code/mywork

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/leshui/code/mywork/src

# Include any dependencies generated for this target.
include CMakeFiles/mywork.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mywork.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mywork.dir/flags.make

CMakeFiles/mywork.dir/main.cpp.o: CMakeFiles/mywork.dir/flags.make
CMakeFiles/mywork.dir/main.cpp.o: main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leshui/code/mywork/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mywork.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mywork.dir/main.cpp.o -c /home/leshui/code/mywork/src/main.cpp

CMakeFiles/mywork.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mywork.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leshui/code/mywork/src/main.cpp > CMakeFiles/mywork.dir/main.cpp.i

CMakeFiles/mywork.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mywork.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leshui/code/mywork/src/main.cpp -o CMakeFiles/mywork.dir/main.cpp.s

CMakeFiles/mywork.dir/material.cpp.o: CMakeFiles/mywork.dir/flags.make
CMakeFiles/mywork.dir/material.cpp.o: material.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leshui/code/mywork/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/mywork.dir/material.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mywork.dir/material.cpp.o -c /home/leshui/code/mywork/src/material.cpp

CMakeFiles/mywork.dir/material.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mywork.dir/material.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leshui/code/mywork/src/material.cpp > CMakeFiles/mywork.dir/material.cpp.i

CMakeFiles/mywork.dir/material.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mywork.dir/material.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leshui/code/mywork/src/material.cpp -o CMakeFiles/mywork.dir/material.cpp.s

CMakeFiles/mywork.dir/vec3.cpp.o: CMakeFiles/mywork.dir/flags.make
CMakeFiles/mywork.dir/vec3.cpp.o: vec3.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leshui/code/mywork/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/mywork.dir/vec3.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mywork.dir/vec3.cpp.o -c /home/leshui/code/mywork/src/vec3.cpp

CMakeFiles/mywork.dir/vec3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mywork.dir/vec3.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leshui/code/mywork/src/vec3.cpp > CMakeFiles/mywork.dir/vec3.cpp.i

CMakeFiles/mywork.dir/vec3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mywork.dir/vec3.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leshui/code/mywork/src/vec3.cpp -o CMakeFiles/mywork.dir/vec3.cpp.s

CMakeFiles/mywork.dir/hittable.cpp.o: CMakeFiles/mywork.dir/flags.make
CMakeFiles/mywork.dir/hittable.cpp.o: hittable.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leshui/code/mywork/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/mywork.dir/hittable.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mywork.dir/hittable.cpp.o -c /home/leshui/code/mywork/src/hittable.cpp

CMakeFiles/mywork.dir/hittable.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mywork.dir/hittable.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leshui/code/mywork/src/hittable.cpp > CMakeFiles/mywork.dir/hittable.cpp.i

CMakeFiles/mywork.dir/hittable.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mywork.dir/hittable.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leshui/code/mywork/src/hittable.cpp -o CMakeFiles/mywork.dir/hittable.cpp.s

CMakeFiles/mywork.dir/sphere.cpp.o: CMakeFiles/mywork.dir/flags.make
CMakeFiles/mywork.dir/sphere.cpp.o: sphere.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leshui/code/mywork/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/mywork.dir/sphere.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mywork.dir/sphere.cpp.o -c /home/leshui/code/mywork/src/sphere.cpp

CMakeFiles/mywork.dir/sphere.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mywork.dir/sphere.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leshui/code/mywork/src/sphere.cpp > CMakeFiles/mywork.dir/sphere.cpp.i

CMakeFiles/mywork.dir/sphere.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mywork.dir/sphere.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leshui/code/mywork/src/sphere.cpp -o CMakeFiles/mywork.dir/sphere.cpp.s

CMakeFiles/mywork.dir/hittable_list.cpp.o: CMakeFiles/mywork.dir/flags.make
CMakeFiles/mywork.dir/hittable_list.cpp.o: hittable_list.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/leshui/code/mywork/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/mywork.dir/hittable_list.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mywork.dir/hittable_list.cpp.o -c /home/leshui/code/mywork/src/hittable_list.cpp

CMakeFiles/mywork.dir/hittable_list.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mywork.dir/hittable_list.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/leshui/code/mywork/src/hittable_list.cpp > CMakeFiles/mywork.dir/hittable_list.cpp.i

CMakeFiles/mywork.dir/hittable_list.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mywork.dir/hittable_list.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/leshui/code/mywork/src/hittable_list.cpp -o CMakeFiles/mywork.dir/hittable_list.cpp.s

# Object files for target mywork
mywork_OBJECTS = \
"CMakeFiles/mywork.dir/main.cpp.o" \
"CMakeFiles/mywork.dir/material.cpp.o" \
"CMakeFiles/mywork.dir/vec3.cpp.o" \
"CMakeFiles/mywork.dir/hittable.cpp.o" \
"CMakeFiles/mywork.dir/sphere.cpp.o" \
"CMakeFiles/mywork.dir/hittable_list.cpp.o"

# External object files for target mywork
mywork_EXTERNAL_OBJECTS =

mywork: CMakeFiles/mywork.dir/main.cpp.o
mywork: CMakeFiles/mywork.dir/material.cpp.o
mywork: CMakeFiles/mywork.dir/vec3.cpp.o
mywork: CMakeFiles/mywork.dir/hittable.cpp.o
mywork: CMakeFiles/mywork.dir/sphere.cpp.o
mywork: CMakeFiles/mywork.dir/hittable_list.cpp.o
mywork: CMakeFiles/mywork.dir/build.make
mywork: CMakeFiles/mywork.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/leshui/code/mywork/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX executable mywork"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mywork.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mywork.dir/build: mywork

.PHONY : CMakeFiles/mywork.dir/build

CMakeFiles/mywork.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mywork.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mywork.dir/clean

CMakeFiles/mywork.dir/depend:
	cd /home/leshui/code/mywork/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/leshui/code/mywork /home/leshui/code/mywork /home/leshui/code/mywork/src /home/leshui/code/mywork/src /home/leshui/code/mywork/src/CMakeFiles/mywork.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mywork.dir/depend

