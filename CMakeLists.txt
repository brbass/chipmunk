cmake_minimum_required (VERSION 2.8)

project(Chipmunk)

# compiler options
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_VERBOSE_MAKEFILE off)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_INSTALL_PREFIX /home/brbass/research/chipmunk/bin)

# standard flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -o3")

file(GLOB chipmunk_src "src/*.cc")

# set(chipmunk_src src/Check.cc src/Chipmunk.cc src/Gauss_Legendre.cc src/Homogeneous_Problem.cc src/Parser.cc src/Sn_Stochastic.cc src/Sn_Transport.cc)

include_directories(src)

add_executable(chipmunk ${chipmunk_src})

install(TARGETS chipmunk DESTINATION .)
