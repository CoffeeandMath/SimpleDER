cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_COMPILER_FORCED TRUE)
set(CMAKE_CXX_COMPILER_FORCED TRUE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
SET(CMAKE_CXX_FLAGS "-std=c++0x")
project(kevin CXX)


 INCLUDE_DIRECTORIES(
  include
  ${CMAKE_SOURCE_DIR}/include
  )

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})

find_package (Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}/unsupported/)
include_directories(${EIGEN3_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}/supported/)

find_package(OpenMP)




ADD_LIBRARY(quaternion ${CMAKE_SOURCE_DIR}/src/quaternion.cpp)
ADD_LIBRARY(node ${CMAKE_SOURCE_DIR}/src/Node.cpp)
ADD_LIBRARY(frame ${CMAKE_SOURCE_DIR}/src/Frame.cpp)
ADD_LIBRARY(element ${CMAKE_SOURCE_DIR}/src/Element.cpp)
ADD_LIBRARY(material ${CMAKE_SOURCE_DIR}/src/Material.cpp)
ADD_LIBRARY(localconstructor ${CMAKE_SOURCE_DIR}/src/LocalConstructor.cpp)
ADD_LIBRARY(constraint ${CMAKE_SOURCE_DIR}/src/Constraint.cpp)
ADD_LIBRARY(inextensibility ${CMAKE_SOURCE_DIR}/src/Inextensibility.cpp)
ADD_LIBRARY(rod ${CMAKE_SOURCE_DIR}/src/Rod.cpp)
ADD_LIBRARY(solver ${CMAKE_SOURCE_DIR}/src/Solver.cpp)
ADD_LIBRARY(read_file ${CMAKE_SOURCE_DIR}/src/read_file.cpp)

ADD_EXECUTABLE(run ${CMAKE_SOURCE_DIR}/main.cpp)
if(OpenMP_CXX_FOUND)
    TARGET_LINK_LIBRARIES(run OpenMP::OpenMP_CXX)
    TARGET_LINK_LIBRARIES(rod OpenMP::OpenMP_CXX)
endif()
TARGET_LINK_LIBRARIES(run quaternion)
TARGET_LINK_LIBRARIES(run node)
TARGET_LINK_LIBRARIES(frame quaternion)
TARGET_LINK_LIBRARIES(frame node)
TARGET_LINK_LIBRARIES(run frame)
TARGET_LINK_LIBRARIES(element quaternion)
TARGET_LINK_LIBRARIES(element node)
TARGET_LINK_LIBRARIES(element frame)
TARGET_LINK_LIBRARIES(run element)
TARGET_LINK_LIBRARIES(material element)
TARGET_LINK_LIBRARIES(run material)
TARGET_LINK_LIBRARIES(localconstructor material)
TARGET_LINK_LIBRARIES(localconstructor element)
TARGET_LINK_LIBRARIES(run localconstructor)
TARGET_LINK_LIBRARIES(rod node)
TARGET_LINK_LIBRARIES(rod frame)
TARGET_LINK_LIBRARIES(rod element)
TARGET_LINK_LIBRARIES(rod material)
TARGET_LINK_LIBRARIES(rod localconstructor)
TARGET_LINK_LIBRARIES(inextensibility constraint)
TARGET_LINK_LIBRARIES(inextensibility frame)
TARGET_LINK_LIBRARIES(solver rod)
TARGET_LINK_LIBRARIES(solver constraint)
TARGET_LINK_LIBRARIES(run rod)
TARGET_LINK_LIBRARIES(run inextensibility)
TARGET_LINK_LIBRARIES(run constraint)
TARGET_LINK_LIBRARIES(run solver)
TARGET_LINK_LIBRARIES(run read_file)
TARGET_LINK_LIBRARIES(run gsl)




# Define a distclean target to remove every generated file:
ADD_CUSTOM_TARGET(distclean
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
  COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
  COMMAND ${CMAKE_COMMAND} -E remove CMakeCache.txt cmake_install.cmake Makefile
  COMMENT "distclean invoked"
  )