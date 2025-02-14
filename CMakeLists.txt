# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project(cgalexec)


# CGAL and its components
find_package(CGAL QUIET COMPONENTS Qt5)

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

# include for local directory

# include for local package


# Creating entries for target: cgalexec
# ############################

add_executable( cgalexec  cgalexec.cpp )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS cgalexec )

# Link the executable to CGAL and third-party libraries
target_link_libraries(cgalexec PRIVATE CGAL::CGAL )

if (CGAL_Qt5_FOUND)
  add_definitions(-DCGAL_USE_BASIC_VIEWER)
  target_link_libraries(cgalexec PRIVATE CGAL::CGAL_Qt5)
endif()