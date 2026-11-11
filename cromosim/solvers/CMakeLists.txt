cmake_minimum_required(VERSION 3.20)

project(optimization VERSION 0.1.0 LANGUAGES CXX)

include(CTest)
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

option(OPTIMIZATION_BUILD_PYTHON "Build the optional Python extension" OFF)
option(OPTIMIZATION_BUILD_BENCHMARKS "Build the reproducible benchmark" ON)
option(OPTIMIZATION_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(OPTIMIZATION_INSTALL_CPP "Install the standalone C++ library and CMake package" ON)
set(OPTIMIZATION_PYTHON_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}" CACHE STRING
    "Install directory for the Python extension")

add_library(optimization_cpp src/optimization/minimizer.cpp)
add_library(optimization::core ALIAS optimization_cpp)
target_compile_features(optimization_cpp PUBLIC cxx_std_17)
set_target_properties(optimization_cpp PROPERTIES
  EXPORT_NAME core
  CXX_EXTENSIONS OFF
  POSITION_INDEPENDENT_CODE ON
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN YES)
target_compile_definitions(optimization_cpp PRIVATE OPTIMIZATION_BUILDING_LIBRARY)
if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(optimization_cpp PUBLIC OPTIMIZATION_STATIC_DEFINE)
endif()
target_include_directories(optimization_cpp
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

find_package(Threads REQUIRED)
target_link_libraries(optimization_cpp PUBLIC Threads::Threads)

if(MSVC)
  target_compile_options(optimization_cpp PRIVATE /W4 /permissive-)
  if(OPTIMIZATION_WARNINGS_AS_ERRORS)
    target_compile_options(optimization_cpp PRIVATE /WX)
  endif()
else()
  target_compile_options(optimization_cpp PRIVATE -Wall -Wextra -Wpedantic)
  if(OPTIMIZATION_WARNINGS_AS_ERRORS)
    target_compile_options(optimization_cpp PRIVATE -Werror)
  endif()
endif()

if(OPTIMIZATION_BUILD_PYTHON)
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(optimization_python MODULE src/optimization/optimization.cpp)
  target_link_libraries(optimization_python PRIVATE optimization::core)
  target_compile_features(optimization_python PRIVATE cxx_std_17)
  set_target_properties(optimization_python PROPERTIES OUTPUT_NAME _optimization CXX_EXTENSIONS OFF)
  install(TARGETS optimization_python
    LIBRARY DESTINATION "${OPTIMIZATION_PYTHON_INSTALL_DIR}"
    RUNTIME DESTINATION "${OPTIMIZATION_PYTHON_INSTALL_DIR}")
endif()

if(BUILD_TESTING)
  add_executable(optimization_tests tests/test_minimizer.cpp)
  target_link_libraries(optimization_tests PRIVATE optimization::core)
  target_compile_features(optimization_tests PRIVATE cxx_std_17)
  set_target_properties(optimization_tests PROPERTIES CXX_EXTENSIONS OFF)
  add_test(NAME optimization.unit COMMAND optimization_tests)
  set_tests_properties(optimization.unit PROPERTIES TIMEOUT 30)
endif()

if(OPTIMIZATION_BUILD_BENCHMARKS)
  add_executable(optimization_benchmark benchmarks/benchmark.cpp)
  target_link_libraries(optimization_benchmark PRIVATE optimization::core)
  target_compile_features(optimization_benchmark PRIVATE cxx_std_17)
  set_target_properties(optimization_benchmark PROPERTIES CXX_EXTENSIONS OFF)
  if(BUILD_TESTING)
    add_test(NAME optimization.benchmark_smoke COMMAND optimization_benchmark
      --particles 512 --iterations 10 --warmup 1 --repetitions 3 --threads 2)
    set_tests_properties(optimization.benchmark_smoke PROPERTIES
      LABELS "performance;smoke" TIMEOUT 60)
  endif()
endif()

if(OPTIMIZATION_INSTALL_CPP)
  install(TARGETS optimization_cpp
    EXPORT optimizationTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(DIRECTORY include/optimization DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  set(OPTIMIZATION_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/optimization")
  configure_package_config_file(cmake/optimizationConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/optimizationConfig.cmake"
    INSTALL_DESTINATION "${OPTIMIZATION_CMAKE_INSTALL_DIR}")
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/optimizationConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)
  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/optimizationConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/optimizationConfigVersion.cmake"
    DESTINATION "${OPTIMIZATION_CMAKE_INSTALL_DIR}")
  install(EXPORT optimizationTargets
    FILE optimizationTargets.cmake
    NAMESPACE optimization::
    DESTINATION "${OPTIMIZATION_CMAKE_INSTALL_DIR}")
endif()
