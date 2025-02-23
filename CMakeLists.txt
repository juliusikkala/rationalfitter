cmake_minimum_required(VERSION 3.19)

project(PolynomialPinner LANGUAGES C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
   add_compile_options(-fdiagnostics-color=always)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   add_compile_options(-fcolor-diagnostics)
endif()

add_executable(polynomialpinner polynomialpinner.cc polynomial.cc rational.cc matrix.cc optimization.cc)
set_property(TARGET polynomialpinner PROPERTY CXX_STANDARD 17)
set_property(TARGET polynomialpinner PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET polynomialpinner PROPERTY CXX_EXTENSIONS OFF)
