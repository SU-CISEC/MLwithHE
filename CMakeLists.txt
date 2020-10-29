cmake_minimum_required (VERSION 3.5.1)

project(MLwithHE CXX)
set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_C_FLAGS_DEBUG "-D_DEBUG")

#if(APPLE)
#    find_package(OpenMP)
#    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
#    set(OpenMP_CXX_LIB_NAMES "omp")
#    set(OpenMP_omp_LIBRARY /usr/local/opt/libomp/lib/libomp.dylib)
#else()
#    find_package(OpenMP)
#    IF(OPENMP_FOUND)
#        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS}")
#        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#    ELSE()
#        message("Missed OpenMP!")
#    ENDIF()
#endif()


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

## PALISADE EXECUTABLES
add_executable(svm-palisade svm-HE-test.cpp)
add_executable(svm-HE svm-HE.cpp)
add_executable( svmOVR svm.cpp)
