# Building emDNA
emDNA supports the following systems:

- macOS 10.13 or later
- Windows 8.1 or later
- Ubuntu 20.04 LTS (any recent Linux distributions should be fine but this is not tested)

The code base has virtually no system-dependent code so emDNA can most likely be built on any platforms that provide a recent C++ compiler (with C++11/14 support) but this is not tested nor validated.


## Table of Contents
- [Requirements](#requirements)
- [Cloning the repository](#cloning-the-repository)
- [Build instructions](#build-instructions)
    * [macOS and Ubuntu (Linux)](#macos-and-ubuntu--linux-)
    * [Windows installation](#windows-installation)
- [Running the tests](#running-the-tests)


## Requirements
- A recent C++ toolchain is required.
- [CMake](https://cmake.org/) v3.18.0 or later.

Note that, emDNA has only been tested on the following platform/compiler setups:

- Apple Clang on macOS
- GCC on Ubuntu
- MSVC (cl.exe) on Windows

It is very likely that the project can be built on other setups (e.g. GCC+MinGW on Windows) but this has not been tested.


## Cloning the repository
emDNA includes various dependencies as submodule so the `--recursive` option is recommended when cloning the repository:

```
git clone --recursive git@github.com:rty10/emDNA.git emDNA.git
```


## Build instructions
The following notes assume the emDNA repository was cloned in the directory `/somewhere/emDNA.git` (`C:\somewhere\emDNA.git` on Windows) so adapt the various paths used in the commands.

Note: because the project is based on CMake it is possible to pass additional compiler flags if desired but this is not covered in the instructions.

### macOS and Ubuntu (Linux)
1. Create a build and installation directory:

    ```
    $ cd /somewhere
    $ mkdir emDNA-build emDNA-install.
    ```

2. Configure the CMake project:

    ```
    $ cmake -S emDNA.git/ -B emDNA-build/ -DCMAKE_BUILD_TYPE=Release
    ```

    You can change the value of `-DCMAKE_BUILD_TYPE` but the recommended one is `Release`.

3. Build the project

    ```
    cmake --build emDNA-build/ --config Release --parallel 2
    ```

    The value 2 used for `--parallel` indicates the number of jobs to use for the build process and it can be increased based on the number of cores/CPUs available on the system.

4. Install:

    ```
    cmake --install emDNA-build/ --config Release --prefix emDNA-install/
    ```

    This will install the various emDNA command-line tools in `/somewhere/emDNA-install/bin/Release` (or whatever build type was selected during configuration). The tools can then be moved/copied elsewhere.

5. Cleanup: the build directory can be removed once the install step was completed successfully.

### Windows
Coming soon ...


## Running the tests
The building and installation process creates two executables that run a set of unit tests on different components of the software:

- `DNASim-tests`: unit tests for the DNASim library which provides various facilities used in the emDNA tools.
- `emDNA-tests`: unit tests for the emDNA implementation itself.

It is recommended to run these two executables and make sure all tests are successful (e.g. no failure messages).
