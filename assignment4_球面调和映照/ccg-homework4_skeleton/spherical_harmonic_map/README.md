# Spherical Harmonic Map

This C++ project framework is used to help students to implement spherical harmonic map algorithm. It contains a simple opengl viewer.

## Dependencies
 
1. `MeshLib`, a mesh library based on halfedge data structure.
2. `freeglut`, a free-software/open-source alternative to the OpenGL Utility Toolkit (GLUT) library.

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

4. Open the \*.sln using Visual Studio, and complie the solution to test everything is OK or not.

5. Finish your code in your IDE.

6. Compile the project `INSTALL` while you finish your code.
> *One need to copy the freeglut.dll or x64/freeglut.dll(depend on your solution platform) into the folder `bin`,
> if your program cannot find the freeglut.dll.*

7. Run the executable program.
> E.x. 
> ``` bash
> cd bin
> ./SphericalHarmonicMap.exe ../data/bimba.m
> ```

7. Press '?' when your mouse is focused on the glut window, and follow the instruction in the command line window.
> If you can see the following results, then it means that you have finished the spherical harmonic map algorithm. 
> 
> ![A face](../resources/bimba.png) ![Harmonic map](../resources/bimba_harmonic_map.png)

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
> ./SphericalHarmonicMap ../data/bimba.m
> ```
