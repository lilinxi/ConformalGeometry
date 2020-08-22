project(OT2d)

if(APPLE)
  add_definitions(-DMAC_OS)
endif(APPLE)

# Include the 'include' of the project itself
include_directories("include")
# Include MeshLib core and Eigen in 3rdparty
include_directories("${MeshLib_DIR}/core"
                    "${CMAKE_SOURCE_DIR}/3rdparty/eigen"
		    "${CMAKE_SOURCE_DIR}/3rdparty/detri2/include")

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

# Add an executable target called MyDemo to be build from 
# the source files.
add_executable(OT2d ${SRCS})

# Link the libraries of freeglut
if(MSVC)
  if(CMAKE_CL_64)
    target_link_libraries(OT2d "${freeglut_DIR}/lib/x64/freeglut.lib"
                         debug "${CMAKE_SOURCE_DIR}/3rdparty/detri2/lib/windows/x64/detri2d.lib"
                     optimized "${CMAKE_SOURCE_DIR}/3rdparty/detri2/lib/windows/x64/detri2.lib")
  else(CMAKE_CL_64)
    target_link_libraries(OT2d "${freeglut_DIR}/lib/freeglut.lib"
                         debug "${CMAKE_SOURCE_DIR}/3rdparty/detri2/lib/windows/Win32/detri2d.lib"
                     optimized "${CMAKE_SOURCE_DIR}/3rdparty/detri2/lib/windows/Win32/detri2.lib")
  endif(CMAKE_CL_64)
else(MSVC)
  target_link_libraries(OT2d "${OPENGL_LIBRARIES}"
                             "${GLUT_LIBRARY}")
  if(APPLE)
    target_link_libraries(OT2d debug "${CMAKE_SOURCE_DIR}/3rdparty/detri2/lib/mac_os/libdetri2d.dylib"
                           optimized "${CMAKE_SOURCE_DIR}/3rdparty/detri2/lib/mac_os/libdetri2.dylib")
  else(APPLE)
    target_link_libraries(OT2d debug "${CMAKE_SOURCE_DIR}/3rdparty/detri2/lib/linux/libdetri2d.so"
                           optimized "${CMAKE_SOURCE_DIR}/3rdparty/detri2/lib/linux/libdetri2.so")
  endif(APPLE)
endif(MSVC)

# Install the executeable program in the bin folder of the
# current project directory.
install(TARGETS OT2d DESTINATION ${CMAKE_SOURCE_DIR}/bin)
