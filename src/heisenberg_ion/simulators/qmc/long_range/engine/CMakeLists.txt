cmake_minimum_required(VERSION 3.31)
project(cpp_qmc)

set(CMAKE_CXX_STANDARD 20)

add_executable(cpp_qmc main.cpp
        ConfigurationGenerator.cpp
        Estimators.cpp
        ProbabilityTables.cpp
        SimulationParameters.cpp
        VertexTypes.cpp)

target_include_directories(cpp_qmc PRIVATE "include/")
