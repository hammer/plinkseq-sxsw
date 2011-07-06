cmake_minimum_required(VERSION 2.8)

include_directories (${PLINKSEQ_SOURCE_DIR}/lib)
link_directories (${PLINKSEQ_BINARY_DIR}/lib) 

include_directories(${PLINKSEQ_SOURCE_DIR}/dependencies/include)
link_directories(${PLINKSEQ_SOURCE_DIR}/dependencies/lib)

add_executable (pseq 
func.h 
summaries.h 
views.h 
assoc.h 
genic.h 
compare.h 
ibd.h 
ibs.h
loaders.h
clusters.h 
netassoc.h
main.cpp 
func.cpp 
summaries.cpp 
views.cpp 
assoc.cpp 
genic.cpp
compare.cpp 
ibd.cpp 
ibs.cpp 
loaders.cpp 
clusters.cpp 
varfunc.cpp
util.cpp 
glm-tests.cpp 
stepup.cpp 
kbac.cpp 
netassoc.cpp 
gseq.cpp
)

target_link_libraries (pseq plinkseq) 
target_link_libraries (pseq protobuf)
target_link_libraries (pseq z)
target_link_libraries (pseq pthread)
target_link_libraries (pseq dl)
