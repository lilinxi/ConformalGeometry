cmake_minimum_required(VERSION 3.6)

# Name the project
project(CCGHomework)

# Find MeshLib
find_package(MeshLib REQUIRED
             HINTS "${CMAKE_SOURCE_DIR}/3rdparty")
if (NOT MeshLib_FOUND)
    message(FATAL_ERROR "MeshLib Not Found!")
endif (NOT MeshLib_FOUND)

# Find freeglut
find_package(freeglut REQUIRED 
             HINTS "${CMAKE_SOURCE_DIR}/3rdparty")
if (NOT freeglut_FOUND)
    message(FATAL_ERROR "freeglut Not Found!")
endif (NOT freeglut_FOUND)

add_subdirectory("cutgraph")

