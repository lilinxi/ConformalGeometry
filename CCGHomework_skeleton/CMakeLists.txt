cmake_minimum_required(VERSION 3.15)
project(CCGHomework)

set(CMAKE_CXX_STANDARD 14)

#list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

# Find MeshLib
find_package(MeshLib REQUIRED HINTS "${CMAKE_SOURCE_DIR}/3rdparty")
if (NOT MeshLib_FOUND)
    message(FATAL_ERROR "MeshLib Not Found!")
endif (NOT MeshLib_FOUND)
include_directories("${MeshLib_DIR}/core")
message(STATUS "MeshLib library status:")
message(STATUS "    config: ${MeshLib_DIR}")
message(STATUS "    include path: ${MeshLib_DIR}/core")

# Find freeglut
# brew cask install xquartz
# brew install freeglut
find_package(FREEGLUT REQUIRED)
if (NOT FREEGLUT_FOUND)
    message(FATAL_ERROR "freeglut Not Found!")
endif (NOT FREEGLUT_FOUND)
message(STATUS "FREEGLUT library status:")
message(STATUS "    config: ${FREEGLUT_DIR}")
message(STATUS "    libraries: ${FREEGLUT_LIBRARIES}")
message(STATUS "    libraries: ${FREEGLUT_LIBRARY_DIRS}")
message(STATUS "    include path: ${FREEGLUT_INCLUDE_DIR}")
include_directories(${FREEGLUT_INCLUDE_DIR})

include_directories("include")
add_subdirectory("cutgraph")

