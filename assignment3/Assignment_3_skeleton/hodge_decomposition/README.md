# Harmonic Map

This C++ project framework is used to help students to implement holomorphic 1-form algorithm. It contains a simple opengl viewer.

## Dependencies
 
1. `MeshLib`, a mesh library based on halfedge data structure.
2. `freeglut`, a free-software/open-source alternative to the OpenGL Utility Toolkit (GLUT) library.
3. `Eigen`, a C++ template library for linear algebra.

## Directory Structure

``` txt
include          -- The header files of harmonic map algorithm.
src              -- The source files of harmonic map algorithm. 
CMakeLists.txt   -- CMake configuration file.
```

## Configuration

### Windows

1. Install [CMake](https://cmake.org/download/).

2. Download the source code of the C++ framework.
> E.x. I create a folder `projects` in `C:/`, then unzip the source code there.

3. Configure and generate the project for Visual Studio.

> ``` bash
> cd CCGHomework
> mkdir build
> cd build
> cmake ..
> ```
> *One can also finish this step using CMake GUI.*

4. Open the \*.sln using Visual Studio, and complie the solution.

5. Finish your code in your IDE.

6. Run the executable program.
> E.x. 
> ``` bash
> cd bin
> ./HodgeDecomposition.exe ../data/boy/boy.m ../data/boy/boy.open.m ../textures/checker_512.bmp
> ./HodgeDecomposition.exe ../data/eight/eight.m ../data/eight/eight.open.m ../textures/checker_512.bmp
> ```

7. Press '?' when your mouse is focused on the glut window, and follow the instruction in the command line window.
> If you can see the following results, then it means that you have finished the harmonic map algorithm. 
> 
> ![Holomorphic 1-form](../resources/boy_holo_1_form.png) ![Holomorphic 1-form](../resources/eight_holo_1_form.png)

### Linux & Mac

1. Build and compile the code.

> ``` bash
> cd CCGHomework
> mkdir build
> cd build
> cmake ..
> make && make install
> ```

2. Run the executable program.

> ``` bash
> cd ../bin/
> ./HodgeDecomposition ../data/eight/eight.m ../data/eight/eight.open.m ../textures/checker_512.bmp
> ```
