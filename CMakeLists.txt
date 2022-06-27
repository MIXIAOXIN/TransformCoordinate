cmake_minimum_required(VERSION 3.5)
project(TransformCoordinate)

set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -luuid")
find_package(Boost REQUIRED COMPONENTS filesystem thread)
include_directories(${BOOST_INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

find_package(PCL REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})

# find_package(Proj REQUIRED)
# include_directories(${PROJECT_SOURCE_DIR})
include_directories(.)

add_executable(TransformCoordinate
        main.cpp
        Trajectory.cpp
        CoordinateSystemTrans.cpp)
