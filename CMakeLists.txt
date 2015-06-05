################################################################
## Directories definition.
################################################################
set(SHARE_DIR "${CMAKE_SOURCE_DIR}/config")    # Directory containing test and configuration files.

set(BINDIR "${CMAKE_BINARY_DIR}/bin" )    # Directory where to store the executable.
set(DOCDIR "${CMAKE_BINARY_DIR}/doc" )    # Directory where to store the API documentation.
set(LIBDIR "${CMAKE_BINARY_DIR}/lib" )    # Directory where to store the project shared library.

##--------------- Do not edit below this line ----------------##

cmake_minimum_required(VERSION 2.8)
project(SHAPE_Opt CXX)

################################################################
## Set up some variables.
################################################################
set(TEST "test")    # Name of the test source file to compile.

set(SHAPE_OPT "shapeopt")    # Name of the shared library.

set(LIB_INSTALLDIR "${CMAKE_INSTALL_PREFIX}/lib/")    # Installation directory for the shared library.
set(INC_INSTALLDIR "${CMAKE_INSTALL_PREFIX}/include/${SHAPE_OPT}")    # Installation directory for the header files.
set(BIN_INSTALLDIR "${CMAKE_INSTALL_PREFIX}/bin/")    # Installation directory for the executable.
set(SHARE_INSTALLDIR "${CMAKE_INSTALL_PREFIX}/share/${SHAPE_OPT}")    # Installation directory for test and configuration files.

set(SRCDIR  "${CMAKE_SOURCE_DIR}/src" )    # Directory where to find the library sources.
set(TESTDIR "${CMAKE_SOURCE_DIR}/test")    # Directory where to find the test sources.

file(GLOB TEST_SRC ${TESTDIR}/${TEST}.cc)
file(GLOB SRCS ${SRCDIR}/*.cc)
file(GLOB HDRS ${SRCDIR}/*.h)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${LIBDIR}")
set(EXECUTABLE_OUTPUT_PATH "${BINDIR}")

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall --pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING    # Default build type: Release.
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()

################################################################
## Dependencies and third-party libraries.
################################################################
include_directories(SYSTEM "${CMAKE_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}/include")

find_package(Boost 1.50 REQUIRED COMPONENTS math_tr1)
include_directories(SYSTEM ${Boost_INCLUDE_DIR})

find_package(LIBMESH REQUIRED)
include_directories(SYSTEM ${LIBMESH_INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

FIND_PACKAGE(MPI REQUIRED)
INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})

################################################################
## Documentation generation settings.
################################################################
set(DOXYFILE_SOURCE_DIR "${CMAKE_SOURCE_DIR}/doc")
set(DOXYFILE_OUTPUT_DIR "${DOCDIR}")
set(DOXYFILE_LATEX "ON")

include(RunDoxygen)

################################################################
## Artistic style.
################################################################
find_package(AStyle)

if(ASTYLE_FOUND)
    add_custom_target(astyle ALL
        COMMAND ${ASTYLE_EXECUTABLE} -q -A1 -s4 -C -S -N -Y -f -p -H -E ${SRCS} ${HDRS} ${TEST_SRC}
        COMMENT "Formatting source codes..."
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endif()

################################################################
## Copy test and configuration files; create documentation dir.
################################################################
file(COPY ${SHARE_DIR} DESTINATION ${CMAKE_BINARY_DIR})
file(MAKE_DIRECTORY ${DOCDIR})

################################################################
## Shared library.
################################################################
add_library(${SHAPE_OPT} SHARED ${SRCS})

set_target_properties(${SHAPE_OPT} PROPERTIES OUTPUT_NAME "${SHAPE_OPT}")    # Shared library filename.
target_link_libraries(${SHAPE_OPT} ${Boost_LIBRARIES})

################################################################
## Target.
################################################################
add_executable(${TEST} ${TEST_SRC} ${SRCS})
#add_executable(${TEST} ${TEST_SRC})

set_target_properties(${TEST} PROPERTIES OUTPUT_NAME "${TEST}"    # Executable filename.
                              INSTALL_RPATH "${LIB_INSTALLDIR}")    # rpath after installation.
target_link_libraries(${TEST} ${LIBMESH_LIBRARIES} ${MPI_LIBRARIES} pthread)
#target_link_libraries(${TEST} ${SHAPE_OPT} ${LIBMESH_LIBRARIES} ${MPI_LIBRARIES} pthread)

################################################################
## Installation.
################################################################
install(TARGETS ${SHAPE_OPT} DESTINATION "${LIB_INSTALLDIR}")
install(FILES ${HDRS} DESTINATION "${INC_INSTALLDIR}")

install(DIRECTORY ${SHARE_DIR} DESTINATION "${SHARE_INSTALLDIR}")

install(TARGETS ${TEST} DESTINATION "${BIN_INSTALLDIR}")

################################################################
## Uninstallation.
################################################################
add_custom_target(uninstall
    "${CMAKE_COMMAND}" -P "${CMAKE_MODULE_PATH}/Uninstall.cmake"
)
