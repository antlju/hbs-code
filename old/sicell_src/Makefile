CC=g++
CFLAGS=-c -Wall -std=c++14
LDFLAGS=
HEADERS=Makefile common.h field.h utils.h diffops.h
DIFFSRC=derivatives.cpp curl.cpp div.cpp grad.cpp lapl.cpp
SOURCES=divtest_v1.cpp $(DIFFSRC)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(HEADERS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o: $(HEADERS)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
