cmake_minimum_required(VERSION 3.9)
project(msfitslib_core VERSION 1.0.1 DESCRIPTION "msfitslib_core description")
include(GNUInstallDirs)
add_library(msfitslib_core SHARED 
     bg_bedlam.cpp  bg_fits.cpp  bg_globals.cpp  bg_stat.cpp         bg_units.cpp  calsol_values.cpp weather_station.cpp
     bg_array.cpp  bg_date.cpp    bg_geo.cpp   bg_norm.cpp     bg_total_power.cpp  bg_vis.cpp    cvalue_vector.cpp  libnova_interface.cpp  
     spectrometer.cpp)

add_definitions(-D_UNIX -D_ENABLE_PROFILER_ -D_ENABLE_TRACING_ -D_ENABLE_VISUALIZATION_ -D_MONTE_CARLO_ -D_NO_ROOT_)
set_target_properties(msfitslib_core PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    PUBLIC_HEADER   "msfitslib.h;spectrometer.h;bg_date.h;bg_globals.h;cvalue_vector.h;calsol_values.h;myparser.h;mystring.h;mystrtable.h;mytypes.h;basedefines.h;basestring.h;myfile.h;basestructs.h")
configure_file(msfitslib_core.pc.in msfitslib_core.pc @ONLY)
target_include_directories(msfitslib_core PRIVATE .)
install(TARGETS msfitslib_core
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(FILES ${CMAKE_BINARY_DIR}/msfitslib_core.pc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)


# EXACUTABLES :
find_library(CFITSIO_LIB cfitsio)
find_library(LIBNOVA_LIB nova)
find_path(FITSIO_INCLUDE_DIR NAMES fitsio.h)
find_library(FFTW3_LIB fftw3 REQUIRED)
# set(CMAKE_LIBRARY_PATH ../..)
# find_library(MSCOMMONLIB mscommonlib REQUIRED)
# set(MSCOMMONLIB ../mscommonlib/lib/libmscommonlib.so)
include_directories(${FITSIO_INCLUDE_DIR})

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
# include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(../mscommonlib/src/)


add_executable(nan_test nan_test.cpp)
add_executable(libtest  libtest.cpp)
target_link_libraries(libtest ${CMAKE_CURRENT_SOURCE_DIR}/../mscommonlib/lib/libmscommonlib.so msfitslib_core ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB})
add_executable(radec2azh radec2azh.cpp)
target_link_libraries(radec2azh ${CMAKE_CURRENT_SOURCE_DIR}/../mscommonlib/lib/libmscommonlib.so msfitslib_core ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB})
add_executable(sid2ux sid2ux.cpp)
target_link_libraries(sid2ux ${CMAKE_CURRENT_SOURCE_DIR}/../mscommonlib/lib/libmscommonlib.so msfitslib_core ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB})
add_executable(ux2sid   ux2sid.cpp) 
target_link_libraries(ux2sid ${CMAKE_CURRENT_SOURCE_DIR}/../mscommonlib/lib/libmscommonlib.so msfitslib_core ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB})
add_executable(ux2sid_file   ux2sid_file.cpp)
target_link_libraries(ux2sid_file ${CMAKE_CURRENT_SOURCE_DIR}/../mscommonlib/lib/libmscommonlib.so msfitslib_core ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB})
