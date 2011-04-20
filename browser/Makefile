
CXX = g++
CXXFLAGS = -O2 

LFLAGS += -Wl,-rpath=/psych/genetics/shaun/pseq/lib 
LFLAGS += -Wl,-rpath=/psych/genetics/pseq/share/lib
LFLAGS += -Wl,-rpath=/psych/genetics/pseq/share/lib64/R/lib

OUTPUT = pbrowse.cgi

LIBDIR = ../lib
LIB = ${LIBDIR}/libpseq.so

SRC = psb.cpp 
HDR = cgi.h psb.h
OBJ = $(SRC:.cpp=.o)

INC = -I../lib/
INC += -I../share/include/ 
INC += -I../share/lib64/R/include/ 

LIBDIR = -L../share/lib -L../lib
LIB  = -lpseq 
LIB += -lsqlite3
LIB += -lprotobuf
LIB += -lRmath
LIB += -lz
LIB += -lpthread

all : $(OUTPUT)

$(OUTPUT) :
	$(CXX) $(CXXFLAGS) $(LFLAGS) -o $(OUTPUT) $(OBJ) $(LIBDIR) $(LIB)

$(OBJ) : $(HDR) 

.cpp.o : 
	$(CXX) $(CXXFLAGS) -c $*.cpp $(INC)

.SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)

FORCE:
