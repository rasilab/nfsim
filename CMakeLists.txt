cmake_minimum_required(VERSION 3.16)

project(NFsim)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(SUB_DIRS
    src/nauty24 
    src/NFutil 
    src/NFutil/MTrand 
    src/NFtest/transcription 
    src/NFtest/tlbr 
    src/NFtest/simple_system 
    src/NFtest/agentcell/cell 
    src/NFtest/agentcell 
    src/NFscheduler 
    src/NFreactions/transformations 
    src/NFreactions/reactions 
    src/NFreactions/reactantLists 
    src/NFreactions/mappings 
    src/NFreactions 
    src/NFoutput 
    src/NFinput 
    src/NFinput/TinyXML 
    src/NFinput/json 
    src/NFfunction/muParser 
    src/NFfunction 
    src/NFcore 
    src/NFcore/reactionSelector 
    src/NFcore/moleculeLists 
)

add_compile_options(-Wno-deprecated-declarations -Wno-format)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} include ${SUB_DIRS})


file(GLOB_RECURSE src_files  "${CMAKE_CURRENT_SOURCE_DIR}/src/*cpp")
file(GLOB_RECURSE c_src_files  "${CMAKE_CURRENT_SOURCE_DIR}/src/*c")


set(SRC_FILES 
    ${src_files} 
    ${c_src_files} 
)


add_executable(${PROJECT_NAME} ${SRC_FILES} )
