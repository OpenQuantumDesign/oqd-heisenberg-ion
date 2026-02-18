cmake_minimum_required(VERSION 3.31)
project(cpp_qmc)

set(CMAKE_CXX_STANDARD 20)

include(FetchContent)

FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG v1.17.0
)

set(SPDLOG_BUILD_EXAMPLE OFF)
set(SPDLOG_BUILD_TESTS OFF)
set(SPDLOG_BUILD_BENCH OFF)
FetchContent_MakeAvailable(spdlog)

add_executable(cpp_qmc main.cpp
        ConfigurationGenerator.cpp
        Estimators.cpp
        ProbabilityTables.cpp
        SimulationParameters.cpp
        VertexTypes.cpp)

target_link_libraries(cpp_qmc PRIVATE spdlog::spdlog)
