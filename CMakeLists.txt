cmake_minimum_required(VERSION 4.0)
project(SI2530PROJECT)

set(CMAKE_CXX_STANDARD 20)

add_executable(SI2530PROJECT main.cpp
        qpottsmodel.cpp
        qpottsmodel.h
        wanglandausimulator.h
)
