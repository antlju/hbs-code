CC=g++
CFLAGS=-c -Wall -std=c++14
LDFLAGS=
HEADERS=Makefile common.h field.h vartypedefs.h utils.h
SOURCES=abcfield_v1.cpp pbc.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(HEADERS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o: $(HEADERS)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
