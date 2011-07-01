cmake_minimum_required (VERSION 2.8) 
project (PLINKSEQ) 

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PLINKSEQ_SOURCE_DIR}/cmake/Modules/")

add_subdirectory (dependencies)
add_subdirectory (lib) 
add_subdirectory (client)
#add_subdirectory (acdb)
