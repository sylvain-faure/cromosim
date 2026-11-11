cmake_minimum_required(VERSION 3.20)

project(cromosim VERSION 3.1.2 LANGUAGES CXX)

# Package builds only need the Python extension. The standalone C++ project in
# cromosim/solvers keeps its own tests, benchmark and installable CMake package.
option(BUILD_TESTING "Build native tests" OFF)
option(OPTIMIZATION_BUILD_BENCHMARKS "Build the native benchmark" OFF)
option(OPTIMIZATION_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
set(OPTIMIZATION_BUILD_PYTHON ON CACHE BOOL "Build the Python extension" FORCE)
set(OPTIMIZATION_INSTALL_CPP OFF CACHE BOOL "Install the standalone C++ package" FORCE)
set(OPTIMIZATION_PYTHON_INSTALL_DIR "cromosim" CACHE STRING
    "Wheel-relative directory for the Python extension" FORCE)

add_subdirectory(cromosim/solvers)
