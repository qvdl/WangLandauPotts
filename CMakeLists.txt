cmake_minimum_required(VERSION 4.0)
project(SI2530PROJECT)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenMP REQUIRED)

add_executable(SI2530PROJECT main.cpp
        qpottsmodel.cpp
        qpottsmodel.h
        wanglandausimulator.h
        rng.h
        wanglandausimulator.cpp
        ThermoState.cpp
        ThermoState.h
        functions.cpp
        functions.h
)

target_link_libraries(SI2530PROJECT PRIVATE OpenMP::OpenMP_CXX)