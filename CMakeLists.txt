cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ilp_example)

find_package(PCL 1.3 REQUIRED COMPONENTS common io)
find_package(CGAL REQUIRED COMPONENTS Core)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# add_compile_options(-fopenmp)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)

include_directories(${PCL_INCLUDE_DIRS})
include_directories(/opt/gurobi910/linux64/include/ include/)

link_directories(${PCL_LIBRARY_DIRS} /opt/gurobi910/linux64/lib)

add_definitions(-DCGAL_USE_BASIC_VIEWER)
add_definitions(${PCL_DEFINITIONS})

add_executable(3d_coverage max_quality.cpp cumulative_quality.cpp bare_visibility.cpp Solution.cpp)

target_link_libraries(3d_coverage ${PCL_LIBRARIES} gurobi_c++ gurobi91 m CGAL::CGAL CGAL::CGAL_Core OpenMP::OpenMP_CXX)
