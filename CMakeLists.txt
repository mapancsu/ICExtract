cmake_minimum_required(VERSION 2.8.0)

project(ICExtract)

# set module path for eigen, numpy
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}  "../thirdparty/eigen")
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}  "../thirdparty/openms")

# find eigen library
find_package( Eigen3 )
if (EIGEN3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR} )
else()
  message(FATAL_ERROR "Could not find eigen, which is a required depedency")
endif()

#find_library(OpenMSL OpenMS)
#if (OpenMSL_Found)
#  include_directories(${OpenMS_INCLUDE_DIR} )
#else()
#  message(FATAL_ERROR "Could not find openms, which is a required depedency")
#endif()

find_package(OpenMS)
if (OpenMS_FOUND)
  include_directories(${OpenMS_INCLUDE_DIRECTORIES} )
else()
  message(FATAL_ERROR "Could not find openms, which is a required depedency")
endif()

find_package(PythonLibs REQUIRED)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


if(MSVC)
	add_definitions(-DUNICODE -D_UNICODE -D_CRT_SECURE_NO_DEPRECATE /wd4251 /wd4819 /wd4244)
else()
	SET(CMAKE_INSTALL_RPATH "$ORIGIN/")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -mavx2")
endif()

include_directories(cpp)

add_subdirectory(cpp)
add_subdirectory(python)
