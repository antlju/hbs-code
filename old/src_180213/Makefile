CC=g++
CFLAGS=-c -Wall -std=c++14
LDFLAGS=
HEADERS=Makefile common.h vartypedefs.h pmesh.h diffops.h field.h qjkmap.h utils.h includes.h maps.h
DIFFOPS=grad.cpp derivatives.cpp
SOURCES=curltest_v2.cpp $(DIFFOPS) maps.cpp pbc.cpp curl.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=curltest.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(HEADERS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o: $(HEADERS)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o *.exe
