CC=g++
CFLAGS=-c -Wall -std=c++14 -Iheaders
LDFLAGS=

HEADERS=Makefile headers/*.h

DIFFOPS=partials/partials.cpp findiff/findiff.cpp curl/curl.cpp div/div.cpp grad/grad.cpp lapl/lapl.cpp
UTILS=utility/bundlememcpy.cpp utility/linspace.cpp utility/pbc.cpp utility/printfield.cpp
MAINSRC=rkdev.cpp
SOURCES=$(MAINSRC) $(DIFFOPS) $(UTILS)

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(HEADERS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o: $(HEADERS)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o *.exe
