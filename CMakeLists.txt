cmake_minimum_required(VERSION 3.19)

project(PolynomialPinner LANGUAGES C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
   add_compile_options(-fdiagnostics-color=always)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   add_compile_options(-fcolor-diagnostics)
endif()

#add_compile_options(-fprofile-arcs -ftest-coverage)
#add_link_options(-fprofile-arcs -ftest-coverage)

# Core is built as a "library" so that it's easier to link to test programs.
add_library(fit-core OBJECT src/polynomial.cc src/rational.cc src/matrix.cc src/optimization.cc)
add_executable(rationalfitter src/rationalfitter.cc)
target_link_libraries(rationalfitter PUBLIC fit-core)
set_property(TARGET rationalfitter PROPERTY CXX_STANDARD 17)
set_property(TARGET rationalfitter PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET rationalfitter PROPERTY CXX_EXTENSIONS OFF)

include(CTest)
enable_testing()
add_subdirectory(test)
