# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

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

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/sdi2000260/Project01

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/sdi2000260/Project01

# Include any dependencies generated for this target.
include CMakeFiles/cgalexec.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/cgalexec.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/cgalexec.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cgalexec.dir/flags.make

CMakeFiles/cgalexec.dir/cgalexec.cpp.o: CMakeFiles/cgalexec.dir/flags.make
CMakeFiles/cgalexec.dir/cgalexec.cpp.o: cgalexec.cpp
CMakeFiles/cgalexec.dir/cgalexec.cpp.o: CMakeFiles/cgalexec.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/sdi2000260/Project01/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/cgalexec.dir/cgalexec.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/cgalexec.dir/cgalexec.cpp.o -MF CMakeFiles/cgalexec.dir/cgalexec.cpp.o.d -o CMakeFiles/cgalexec.dir/cgalexec.cpp.o -c /home/sdi2000260/Project01/cgalexec.cpp

CMakeFiles/cgalexec.dir/cgalexec.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/cgalexec.dir/cgalexec.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sdi2000260/Project01/cgalexec.cpp > CMakeFiles/cgalexec.dir/cgalexec.cpp.i

CMakeFiles/cgalexec.dir/cgalexec.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/cgalexec.dir/cgalexec.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sdi2000260/Project01/cgalexec.cpp -o CMakeFiles/cgalexec.dir/cgalexec.cpp.s

# Object files for target cgalexec
cgalexec_OBJECTS = \
"CMakeFiles/cgalexec.dir/cgalexec.cpp.o"

# External object files for target cgalexec
cgalexec_EXTERNAL_OBJECTS =

cgalexec: CMakeFiles/cgalexec.dir/cgalexec.cpp.o
cgalexec: CMakeFiles/cgalexec.dir/build.make
cgalexec: libCGAL_Qt5_moc_and_resources.a
cgalexec: /usr/lib/x86_64-linux-gnu/libgmpxx.so
cgalexec: /usr/lib/x86_64-linux-gnu/libmpfr.so
cgalexec: /usr/lib/x86_64-linux-gnu/libgmp.so
cgalexec: /usr/lib/x86_64-linux-gnu/libQt5OpenGL.so.5.12.8
cgalexec: /usr/lib/x86_64-linux-gnu/libQt5Widgets.so.5.12.8
cgalexec: /usr/lib/x86_64-linux-gnu/libQt5Gui.so.5.12.8
cgalexec: /usr/lib/x86_64-linux-gnu/libQt5Core.so.5.12.8
cgalexec: CMakeFiles/cgalexec.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/sdi2000260/Project01/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable cgalexec"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cgalexec.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cgalexec.dir/build: cgalexec
.PHONY : CMakeFiles/cgalexec.dir/build

CMakeFiles/cgalexec.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cgalexec.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cgalexec.dir/clean

CMakeFiles/cgalexec.dir/depend:
	cd /home/sdi2000260/Project01 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sdi2000260/Project01 /home/sdi2000260/Project01 /home/sdi2000260/Project01 /home/sdi2000260/Project01 /home/sdi2000260/Project01/CMakeFiles/cgalexec.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/cgalexec.dir/depend

