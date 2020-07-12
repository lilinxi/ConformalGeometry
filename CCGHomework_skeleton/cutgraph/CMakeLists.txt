project(CutGraph)

include_directories(include)

#GLOB 会产生一个由所有匹配globbing表达式的文件组成的列表，并将其保存到变量中。
file(GLOB SRCS "include/*.h" "src/*.cpp")

add_executable(CutGraph ${SRCS})

set(FREEGLUT_LIB /usr/local/Cellar/freeglut/3.2.1/lib/libglut.dylib)
target_link_libraries(CutGraph ${FREEGLUT_LIB})

install(TARGETS CutGraph DESTINATION ${CMAKE_SOURCE_DIR}/bin)
