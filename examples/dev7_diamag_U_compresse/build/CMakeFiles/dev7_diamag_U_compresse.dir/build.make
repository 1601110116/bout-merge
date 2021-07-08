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
CMAKE_COMMAND = /vol6/software/cmake-3.14.4/bin/cmake

# The command to remove a file.
RM = /vol6/software/cmake-3.14.4/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /vol6/home/chaodong/bout-nersc/examples/dev7_diamag_U_compresse

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /vol6/home/chaodong/bout-nersc/examples/dev7_diamag_U_compresse/build

# Include any dependencies generated for this target.
include CMakeFiles/dev7_diamag_U_compresse.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/dev7_diamag_U_compresse.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/dev7_diamag_U_compresse.dir/flags.make

CMakeFiles/dev7_diamag_U_compresse.dir/dev7_diamag_U_compresse.cxx.o: CMakeFiles/dev7_diamag_U_compresse.dir/flags.make
CMakeFiles/dev7_diamag_U_compresse.dir/dev7_diamag_U_compresse.cxx.o: ../dev7_diamag_U_compresse.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/vol6/home/chaodong/bout-nersc/examples/dev7_diamag_U_compresse/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/dev7_diamag_U_compresse.dir/dev7_diamag_U_compresse.cxx.o"
	/vol6/software/mpi/mpi-intel2013/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dev7_diamag_U_compresse.dir/dev7_diamag_U_compresse.cxx.o -c /vol6/home/chaodong/bout-nersc/examples/dev7_diamag_U_compresse/dev7_diamag_U_compresse.cxx

CMakeFiles/dev7_diamag_U_compresse.dir/dev7_diamag_U_compresse.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dev7_diamag_U_compresse.dir/dev7_diamag_U_compresse.cxx.i"
	/vol6/software/mpi/mpi-intel2013/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /vol6/home/chaodong/bout-nersc/examples/dev7_diamag_U_compresse/dev7_diamag_U_compresse.cxx > CMakeFiles/dev7_diamag_U_compresse.dir/dev7_diamag_U_compresse.cxx.i

CMakeFiles/dev7_diamag_U_compresse.dir/dev7_diamag_U_compresse.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dev7_diamag_U_compresse.dir/dev7_diamag_U_compresse.cxx.s"
	/vol6/software/mpi/mpi-intel2013/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /vol6/home/chaodong/bout-nersc/examples/dev7_diamag_U_compresse/dev7_diamag_U_compresse.cxx -o CMakeFiles/dev7_diamag_U_compresse.dir/dev7_diamag_U_compresse.cxx.s

# Object files for target dev7_diamag_U_compresse
dev7_diamag_U_compresse_OBJECTS = \
"CMakeFiles/dev7_diamag_U_compresse.dir/dev7_diamag_U_compresse.cxx.o"

# External object files for target dev7_diamag_U_compresse
dev7_diamag_U_compresse_EXTERNAL_OBJECTS =

dev7_diamag_U_compresse: CMakeFiles/dev7_diamag_U_compresse.dir/dev7_diamag_U_compresse.cxx.o
dev7_diamag_U_compresse: CMakeFiles/dev7_diamag_U_compresse.dir/build.make
dev7_diamag_U_compresse: /vol6/home/chaodong/local/bout-nersc/lib64/libbout++.a
dev7_diamag_U_compresse: /vol6/home/chaodong/local/bout-nersc/lib64/libpvode.a
dev7_diamag_U_compresse: /vol6/home/chaodong/local/bout-nersc/lib64/libpvpre.a
dev7_diamag_U_compresse: /vol6/software/io_tools/netcdf/mpi/4.1.2/lib/libnetcdf_c++4.so
dev7_diamag_U_compresse: /vol6/software/io_tools/netcdf/mpi/4.1.2/lib/libnetcdf.so
dev7_diamag_U_compresse: /vol6/software/libraries/fftw/mpi/3.3.3/lib/libfftw3.so
dev7_diamag_U_compresse: /vol6/home/chaodong/software/petsc-3.5.4/arch-linux2-cxx-debug/lib/libpetsc.so
dev7_diamag_U_compresse: CMakeFiles/dev7_diamag_U_compresse.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/vol6/home/chaodong/bout-nersc/examples/dev7_diamag_U_compresse/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable dev7_diamag_U_compresse"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dev7_diamag_U_compresse.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/dev7_diamag_U_compresse.dir/build: dev7_diamag_U_compresse

.PHONY : CMakeFiles/dev7_diamag_U_compresse.dir/build

CMakeFiles/dev7_diamag_U_compresse.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/dev7_diamag_U_compresse.dir/cmake_clean.cmake
.PHONY : CMakeFiles/dev7_diamag_U_compresse.dir/clean

CMakeFiles/dev7_diamag_U_compresse.dir/depend:
	cd /vol6/home/chaodong/bout-nersc/examples/dev7_diamag_U_compresse/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /vol6/home/chaodong/bout-nersc/examples/dev7_diamag_U_compresse /vol6/home/chaodong/bout-nersc/examples/dev7_diamag_U_compresse /vol6/home/chaodong/bout-nersc/examples/dev7_diamag_U_compresse/build /vol6/home/chaodong/bout-nersc/examples/dev7_diamag_U_compresse/build /vol6/home/chaodong/bout-nersc/examples/dev7_diamag_U_compresse/build/CMakeFiles/dev7_diamag_U_compresse.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/dev7_diamag_U_compresse.dir/depend
