project(SphericalHarmonicMap)

if(APPLE)
  add_definitions(-DMAC_OS)
endif(APPLE)

# Include the 'include' of the project itself
include_directories("include")
# Include MeshLib
include_directories("${MeshLib_DIR}/core")

# Include the directory of OpenGL and GLUT
if(MSVC)
  include_directories("${freeglut_DIR}/include")
else(MSVC)
  include_directories(${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR})
endif(MSVC)

# The file(GLOB...) allows for wildcard additions
file(GLOB SRCS
    "include/*.h"
    "src/*.cpp")

# Add an executable target called SphericalHarmonicMap to be build from 
# the source files.
add_executable(SphericalHarmonicMap ${SRCS})

# Link the libraries of freeglut
if(MSVC)
  if(CMAKE_CL_64)
    target_link_libraries(SphericalHarmonicMap "${freeglut_DIR}/lib/x64/freeglut.lib")
  else(CMAKE_CL_64)
    target_link_libraries(SphericalHarmonicMap "${freeglut_DIR}/lib/freeglut.lib")
  endif(CMAKE_CL_64)
else(MSVC)
  target_link_libraries(SphericalHarmonicMap ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
endif(MSVC)

# Install the executeable program in the bin folder of the
# current project directory.
install(TARGETS SphericalHarmonicMap DESTINATION ${CMAKE_SOURCE_DIR}/bin)
