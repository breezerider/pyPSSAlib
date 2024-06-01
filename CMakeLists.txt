cmake_minimum_required(VERSION 3.15)
project(pypssalib VERSION ${PYPSSALIB_VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(GNUInstallDirs)
include(ExternalProject)

#set(Boost_DEBUG ON)
#set(Boost_NO_BOOST_CMAKE    ON)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED CONFIG)
find_package(GSL REQUIRED)
find_package(Boost 1.56 REQUIRED)

ExternalProject_add(
    pssalib
    BUILD_BYPRODUCTS ${CMAKE_CURRENT_SOURCE_DIR}/externals/${CMAKE_INSTALL_LIBDIR}/libpssa${CMAKE_STATIC_LIBRARY_SUFFIX}
    INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/externals"
    GIT_REPOSITORY https://github.com/breezerider/pSSAlib.git
    GIT_TAG v2.0.1
    PREFIX pssalib
    CONFIGURE_COMMAND libtoolize && aclocal && automake --add-missing && autoreconf --force && ./configure --prefix=${CMAKE_CURRENT_SOURCE_DIR}/externals --libdir=${CMAKE_CURRENT_SOURCE_DIR}/externals/${CMAKE_INSTALL_LIBDIR} --enable-static --disable-shared --with-pic --without-libsbml --disable-cli --disable-optimisation
    BUILD_COMMAND make -j
    INSTALL_COMMAND make install
    BUILD_IN_SOURCE ON
    )
add_library(libpssa STATIC IMPORTED)
set_target_properties(libpssa PROPERTIES
    IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/externals/${CMAKE_INSTALL_LIBDIR}/libpssa${CMAKE_STATIC_LIBRARY_SUFFIX}
)
add_dependencies(libpssa pssalib)

add_library(pypssalib MODULE src/main.cpp)
set_target_properties(pypssalib PROPERTIES
                                INTERPROCEDURAL_OPTIMIZATION ON
                                CXX_VISIBILITY_PRESET default
                                VISIBILITY_INLINES_HIDDEN ON
                                PREFIX ""
                                SUFFIX ".${PYTHON_SOABI}.so")
target_compile_definitions(pypssalib PRIVATE PYPSSALIB_VERSION=${PYPSSALIB_VERSION})
target_include_directories(pypssalib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/externals/include)
target_include_directories(pypssalib PRIVATE ${Python3_INCLUDE_DIRS})
target_include_directories(pypssalib PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(pypssalib pybind11::module)
target_link_libraries(pypssalib libpssa)
target_link_libraries(pypssalib Python3::Module)
target_link_libraries(pypssalib GSL::gsl)
add_dependencies(pypssalib libpssa)
