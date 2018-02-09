CC=g++
CFLAGS=-c -Wall -std=c++14
LDFLAGS=
HEADERS=Makefile common.h vartypedefs.h pmesh.h diffops.h field.h qjkmap.h includes.h utils.h maps.h
DIFFOPS=grad.cpp derivatives.cpp
SOURCES=maptest_v2.cpp $(DIFFOPS)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(HEADERS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o: $(HEADERS)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
