cmake_minimum_required(VERSION 3.20)
project(cme_detect VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(CME_USE_OPENMP "Enable OpenMP parallelism over frames" ON)

# ---------------------------------------------------------------------------
# Dependencies: GSL (distributions, vectors/matrices, BLAS) and libtiff.
# Pre-built binaries for x64-windows/MSVC are vendored under third_party/;
# on other platforms (or with -DCME_USE_SYSTEM_LIBS=ON) find_package is used.
# The Levenberg-Marquardt solver itself is compiled from the GSL 1.16 sources
# in third_party/gsl116 (see PORTING_NOTES.md).
# ---------------------------------------------------------------------------
option(CME_USE_SYSTEM_LIBS "Use find_package(GSL)/find_package(TIFF) instead of third_party/" OFF)
set(TP ${CMAKE_CURRENT_SOURCE_DIR}/third_party)

if(NOT CME_USE_SYSTEM_LIBS AND EXISTS ${TP}/gsl/lib/gsl.lib AND EXISTS ${TP}/tiff/lib/tiff.lib)
  add_library(gsl_dep SHARED IMPORTED)
  set_target_properties(gsl_dep PROPERTIES
    IMPORTED_LOCATION ${TP}/gsl/bin/gsl.dll
    IMPORTED_IMPLIB ${TP}/gsl/lib/gsl.lib
    INTERFACE_INCLUDE_DIRECTORIES ${TP}/gsl/include)
  add_library(gslcblas_dep SHARED IMPORTED)
  set_target_properties(gslcblas_dep PROPERTIES
    IMPORTED_LOCATION ${TP}/gsl/bin/gslcblas.dll
    IMPORTED_IMPLIB ${TP}/gsl/lib/gslcblas.lib)
  add_library(tiff_dep STATIC IMPORTED)
  set_target_properties(tiff_dep PROPERTIES
    IMPORTED_LOCATION ${TP}/tiff/lib/tiff.lib
    INTERFACE_INCLUDE_DIRECTORIES ${TP}/tiff/include)
  set(CME_DEP_LIBS gsl_dep gslcblas_dep tiff_dep)
  set(CME_VENDORED_GSL ON)
else()
  find_package(GSL REQUIRED)
  find_package(TIFF REQUIRED)
  set(CME_DEP_LIBS GSL::gsl GSL::gslcblas TIFF::TIFF)
  set(CME_VENDORED_GSL OFF)
endif()

# GSL 1.16 Levenberg-Marquardt (lmsder) + covariance, compiled as C
add_library(gsl116_lm STATIC
  ${TP}/gsl116/g116_lm.c
  ${TP}/gsl116/g116_linalg.c)
target_include_directories(gsl116_lm PUBLIC ${TP}/gsl116)
target_link_libraries(gsl116_lm PUBLIC ${CME_DEP_LIBS})
# inline GSL vector/matrix accessors (no change in arithmetic, ~2x faster fits)
target_compile_definitions(gsl116_lm PUBLIC HAVE_INLINE GSL_RANGE_CHECK_OFF)
if(MSVC)
  target_compile_definitions(gsl116_lm PRIVATE _CRT_SECURE_NO_WARNINGS)
  target_compile_options(gsl116_lm PRIVATE /fp:precise /wd4244 /wd4267)
else()
  target_compile_options(gsl116_lm PRIVATE -ffp-contract=off)
endif()

add_library(cme_core STATIC
  src/conv.cpp
  src/morphology.cpp
  src/stats.cpp
  src/tiff_io.cpp
  src/fit_gaussian2d.cpp
  src/fit_gaussians2d.cpp
  src/point_source_detection.cpp
  src/gmm.cpp
  src/psf_sigma.cpp
  src/condition_data.cpp
  src/run_detection.cpp
  src/dump.cpp
  src/output.cpp)
target_include_directories(cme_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cme_core PUBLIC ${CME_DEP_LIBS} gsl116_lm)
target_compile_definitions(cme_core PUBLIC _USE_MATH_DEFINES NOMINMAX WIN32_LEAN_AND_MEAN _CRT_SECURE_NO_WARNINGS)
if(MSVC)
  target_compile_options(cme_core PUBLIC /W3 /fp:precise /EHsc)
else()
  target_compile_options(cme_core PUBLIC -Wall -Wextra -ffp-contract=off)
endif()

if(CME_USE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(cme_core PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()

add_executable(cme_detect src/main.cpp)
target_link_libraries(cme_detect PRIVATE cme_core)

add_executable(cme_tests tests/unit_tests.cpp)
target_link_libraries(cme_tests PRIVATE cme_core)

add_executable(fit_probe tests/fit_probe.cpp)
target_link_libraries(fit_probe PRIVATE cme_core)

add_executable(profile_sigma tests/profile_sigma.cpp)
target_link_libraries(profile_sigma PRIVATE cme_core)

if(CME_VENDORED_GSL)
  foreach(t cme_detect cme_tests fit_probe profile_sigma)
    add_custom_command(TARGET ${t} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${TP}/gsl/bin/gsl.dll ${TP}/gsl/bin/gslcblas.dll $<TARGET_FILE_DIR:${t}>)
  endforeach()
endif()

enable_testing()
add_test(NAME unit_tests COMMAND cme_tests ${CMAKE_CURRENT_SOURCE_DIR}/tests/reference_values.txt)
