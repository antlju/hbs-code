CC=g++
CFLAGS=-c -Wall -std=c++14
LDFLAGS=
HEADERS=Makefile common.h vartypedefs.h pmesh.h diffops.h field.h
DIFFOPS=grad.cpp derivatives.cpp
SOURCES=diffoptest_v2.cpp $(DIFFOPS)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(HEADERS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o: $(HEADERS)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
