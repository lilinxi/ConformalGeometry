project(CutGraph)

# Include MeshLib core and the 'include' of the project itsself
include_directories(
    "include"
    "${MeshLib_DIR}/core"
    "${freeglut_DIR}/include")

# Manually add the sources
# set(SOURCES_ONE ${MeshLib_DIR}/core/bmp/RgbImage.cpp)

# Also, the file(GLOB...) allows for wildcard additions
file(GLOB SRCS
    "include/*.h"
    "src/*.cpp")

# Add an executable target called MyDemo to be build from 
# the source files.
add_executable(CutGraph ${SRCS})

# Link the libraries of freeglut
if(CMAKE_CL_64)
  target_link_libraries(CutGraph
    "${freeglut_DIR}/lib/x64/freeglut.lib")
else(CMAKE_CL_64)
  target_link_libraries(CutGraph
    "${freeglut_DIR}/lib/freeglut.lib")
endif(CMAKE_CL_64)

# Install the executeable program in the bin folder of the
# current project directory.
install(TARGETS CutGraph DESTINATION ${CMAKE_SOURCE_DIR}/bin)
