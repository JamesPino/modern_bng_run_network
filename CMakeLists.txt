cmake_minimum_required(VERSION 3.10.2)

set(CMAKE_CXX_STANDARD 17)

project(Network3)

add_definitions(-DRUN_NETWORK_VERSION=\"3.0\")

set(CMAKELISTS_DIR $(CURDIR))
set(INCDIR include )
set(LIBDIR lib)
set(BNG_BINDIR ../bin)

set(CVODE cvode-2.6.0)
set(LIBSOURCE ../libsource )


set(SUB_DIRS 
    src 
    src/model 
    src/pla 
    src/ssa 
    src/util 
    src/model/rateExpressions 
    src/model/reactions 
    src/pla/base 
    src/pla/eRungeKutta 
    src/pla/fEuler 
    src/pla/util 
    src/pla/eRungeKutta/base
    src/pla/eRungeKutta/util 
    src/util/MTrand 
    src/util/rand2
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} include include/cvode include/nvector include/sundials ${SUB_DIRS})

aux_source_directory(src src_files)
aux_source_directory(src/model model_files)
aux_source_directory(src/pla pla_files)
aux_source_directory(src/ssa ssa_files)
aux_source_directory(src/util src_util_files)
aux_source_directory(src/model/rateExpressions rexp_files)
aux_source_directory(src/model/reactions rxn_files)
aux_source_directory(src/pla/base pla_base_files)
aux_source_directory(src/pla/eRungeKutta eRK_files)
aux_source_directory(src/pla/fEuler fEuler_files)
aux_source_directory(src/pla/util pla_util_files)
aux_source_directory(src/pla/eRungeKutta/base eRK_base_files)
aux_source_directory(src/pla/eRungeKutta/util eRK_util_files)
aux_source_directory(src/util/MTrand mtrand_files)
aux_source_directory(src/util/rand2 rand2_files)

set(SRC_FILES 
    ${src_files} 
    ${model_files} 
    ${pla_files} 
    ${ssa_files}  
    ${src_util_files} 
    ${rexp_files} 
    ${rxn_files} 
    ${pla_base_files} 
    ${eRK_files} 
    ${fEuler_files} 
    ${pla_util_files} 
    ${eRK_base_files} 
    ${eRK_util_files} 
    ${mtrand_files} 
    ${rand2_files}
)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.45.0 COMPONENTS program_options)

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
link_libraries(libmathutils.a libmuparser.a libsundials_cvode.a libsundials_nvecserial.a)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(parse_input ${SRC_FILES})
    target_link_libraries(parse_input ${Boost_LIBRARIES})
endif()

#add_custom_target(
#        myCustomTarget
#        COMMENT "Unpacking abc.zip"
#        COMMAND rm -rf ${CVODE}
#        COMMAND pwd
#        COMMAND mkdir -p ${LIBDIR} ${INCDIR} ${BNG_BINDIR}
#        COMMAND tar -xzf ${LIBSOURCE}/${CVODE}.tar.gz
#        WORKING_DIRECTORY ${CMAKELISTS_DIR}/${CVODE}
#        COMMAND ./configure --prefix=${CURDIR} --disable-shared
#        COMMAND  make
#        COMMAND  make install
#        WORKING_DIRECTORY ${CMAKELISTS_DIR}
#        COMMENT "Unpacking abc.zip"
#        VERBATIM)
#add_dependencies(run_network myCustomTarget)