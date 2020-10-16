cmake_minimum_required (VERSION 3.5.1)

project(MLwithHE CXX)
set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_C_FLAGS_DEBUG "-D_DEBUG")
#find_package(SEAL 3.5.8 EXACT REQUIRED)

if(APPLE)
    find_package(OpenMP)
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY /usr/local/opt/libomp/lib/libomp.dylib)
else()
    find_package(OpenMP)
    IF(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    ELSE()
        message("Missed OpenMP!")
    ENDIF()
endif()
file(GLOB SOURCES
        header-folder/*.h
        source-folder/*.cpp
        )


find_package(Palisade)

set( CMAKE_CXX_FLAGS ${PALISADE_CXX_FLAGS} )
set( CMAKE_EXE_LINKER_FLAGS ${PALISADE_EXE_LINKER_FLAGS} )

include_directories( ${OPENMP_INCLUDES} )
include_directories( ${PALISADE_INCLUDE} )
include_directories( ${PALISADE_INCLUDE}/core )
include_directories( ${PALISADE_INCLUDE}/pke )

link_directories( ${PALISADE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )
link_libraries( ${PALISADE_LIBRARIES} )


#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

### ADD YOUR EXECUTABLE(s) HERE
### add_executable( EXECUTABLE-NAME SOURCES )
###
### EXAMPLE:
#add_executable( bm-log_reg bm-log_reg.cpp )
#target_link_libraries(bm-log_reg.cpp SEAL::seal)

## SEAL EXECUTABLES
#add_executable( bm-chi2 benchmark-chi2.cpp)
#target_link_libraries(bm-chi2 SEAL::seal)

## PALISADE EXECUTABLES
add_executable( svm-palisade palisade/svm-palisade.cpp)
add_executable( svm-seal seal/svm-seal.cpp)
add_executable( svmOVR seal/svmWithOVR.cpp)
