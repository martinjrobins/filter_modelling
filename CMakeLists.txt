cmake_minimum_required(VERSION 2.8)
project(AdaptiveDiffusion)

set(Aboria_LOG_LEVEL 1 CACHE STRING "Logging level (1 = least, 3 = most)")
add_definitions(-DABORIA_LOG_LEVEL=${Aboria_LOG_LEVEL})

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/Aboria/cmake"
                    "${CMAKE_CURRENT_SOURCE_DIR}"
                        ${CMAKE_MODULE_PATH})

list(APPEND CMAKE_CXX_FLAGS "-Wno-deprecated -std=c++14")
list(APPEND CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

find_package(Boost 1.50.0 COMPONENTS program_options serialization REQUIRED)
set(filter_libs ${Boost_LIBRARIES})
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)

#find_package(VTK REQUIRED)
#add_definitions(-DHAVE_VTK)
#set(filter_libs ${filter_libs} ${VTK_LIBRARIES} )

option(Aboria_USE_GPERFTOOLS "Use Google Profiling tools" OFF)
if (Aboria_USE_GPERFTOOLS)
    find_package(Gperftools REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
    set(GPERFTOOLS_LIBRARIES "-Wl,--no-as-needed ${GPERFTOOLS_LIBRARIES} -Wl,--as-needed")
    add_definitions(-DHAVE_GPERFTOOLS)
endif()

option(Aboria_USE_OPENMP "Use OpenMP for shared memory parallism" OFF)
if (Aboria_USE_OPENMP)
    find_package(OpenMP REQUIRED)
    add_definitions(-DHAVE_OPENMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


option(Aboria_USE_THRUST "Use CUDA Thrust library" OFF)
if (Aboria_USE_THRUST)
    find_package(Thrust REQUIRED)
    find_package(CUDA REQUIRED)

    option(Aboria_THRUST_USE_ALGORITHMS "Use CUDA Thrust library algorithms" ON)
    option(Aboria_THRUST_TEST_DEVICE_VECTOR "Turn on tests using thrust::device_vector" OFF)
    if (Aboria_THRUST_USE_ALGORITHMS)
        add_definitions(-D__aboria_use_thrust_algorithms__)
    endif()
    #http://stackoverflow.com/questions/34996295/trying-to-get-cuda-7-5-to-work-with-gcc-5-x

    add_definitions(-D_MWAITXINTRIN_H_INCLUDED)
    set(
    CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS};
    --keep
    #--std=c++11
    #-gencode arch=compute_20,code=sm_20
    )
endif()

find_package(Eigen3 REQUIRED)
add_definitions(-DHAVE_EIGEN)

find_package(CGAL REQUIRED)
find_package(GMP REQUIRED)
set(filter_libs ${filter_libs} ${CGAL_LIBRARIES} ${GMP_LIBRARIES})

#export compiler flags for code completion engines
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

set(source_files
    solve_stokes_MAPS.cpp
    filter.cpp
    setup_knots.cpp
    solve_stokes_LMAPS.cpp
    solve_stokes_Compact.cpp
    )

set(filter_include_dirs "Aboria/src" "Aboria/third-party")
set(filter_include_dirs ${filter_include_dirs} ${CXXTEST_INCLUDES} ${EIGEN3_INCLUDE_DIR} ${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} ${CGAL_INCLUDE_DIR} ${GMP_INCLUDE_DIR})

include_directories(${filter_include_dirs})

add_library(filter_modelling SHARED ${source_files})
target_link_libraries(filter_modelling ${filter_libs})

add_executable(check_solutions check_solutions.cpp )
target_link_libraries(check_solutions ${filter_libs})

add_executable(simulation simulation.cpp )
target_link_libraries(simulation filter_modelling )

add_executable(optimc optimc.cpp )
target_link_libraries(optimc filter_modelling )

add_executable(fmm_eval fmm_eval.cpp )
target_link_libraries(fmm_eval filter_modelling )
