# Scratch notes for build+install procedure #

These are scratch notes ... work in progress / need review

The idea is to drive the entire process from CMake in a "modern" way ...

So let's start with assuming we have the following file structure:

```
somewhere/
    emDNA.git/          <-- the repository
```

0. Clone the repository:

   ```
   cd somewhere/
   git clone --recursive git@github.com:rty10/emDNA.git emDNA.git
   ```

1. Create a build directory (stay in `somewhere/` directory):

    ```
    cd somewhere/
    mkdir build
    ```

2. Configure the cmake project (`Release` can be any of `Debug`, `Release`, `RelWithDebInfo` or `MinSizeRel`):

    ```
    cmake -S emDNA.git -B build/ -DCMAKE_BUILD_TYPE=Release
    ```

3. Build the project:

    ```
    cmake --build build/ --config Release --parallel 2 --clean-first
    ```

    The value `parallel` option can be changed to accommodate how many jobs/cores/CPUs will be used to build the project.

    The `clean-first` option is useful when rebuilding in the same build directory.

3. Install:

    ```
    cmake --install build/ --config Release --prefix install/
    ```

    This will install in `somewhere/install/` but it can be anything. If not defined the install directory is (I believe) platform dependent ... `/usr/local` on Linux/macOS, something else on Windows.

