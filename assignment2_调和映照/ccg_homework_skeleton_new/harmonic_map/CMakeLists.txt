project(HarmonicMap)

if(APPLE)
  add_definitions(-DMAC_OS)
endif(APPLE)

# Include the 'include' of the project itself
include_directories("include")
# Include MeshLib core and Eigen in 3rdparty
include_directories("${MeshLib_DIR}/core"
                    "${CMAKE_SOURCE_DIR}/3rdparty/eigen")

# Include the directory of OpenGL and GLUT
if(MSVC)
  include_directories("${freeglut_DIR}/include")
else(MSVC)
  include_directories(${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR})
endif(MSVC)

# Manually add the sources
# set(SOURCES_ONE ${MeshLib_DIR}/core/bmp/RgbImage.cpp)

# Also, the file(GLOB...) allows for wildcard additions
file(GLOB SRCS
    "include/*.h"
    "src/*.cpp")

# Add an executable target called MyDemo to be build from 
# the source files.
add_executable(HarmonicMap ${SRCS})

# Link the libraries of freeglut
if(MSVC)
  if(CMAKE_CL_64)
    target_link_libraries(HarmonicMap "${freeglut_DIR}/lib/x64/freeglut.lib")
  else(CMAKE_CL_64)
    target_link_libraries(HarmonicMap "${freeglut_DIR}/lib/freeglut.lib")
  endif(CMAKE_CL_64)
else(MSVC)
  target_link_libraries(HarmonicMap ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
endif(MSVC)

# Install the executeable program in the bin folder of the
# current project directory.
install(TARGETS HarmonicMap DESTINATION ${CMAKE_SOURCE_DIR}/bin)
