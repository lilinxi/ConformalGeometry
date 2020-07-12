project(CutGraph)

#GLOB 会产生一个由所有匹配globbing表达式的文件组成的列表，并将其保存到变量中。
file(GLOB SRCS "include/*.h" "src/*.cpp")

add_executable(CutGraph ${SRCS})

target_link_libraries(CutGraph FreeGLUT::freeglut)

install(TARGETS CutGraph DESTINATION ${CMAKE_SOURCE_DIR}/bin)
