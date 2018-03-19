CC=g++
CFLAGS=-c -Wall -std=c++14 -Iheaders -pg
LDFLAGS=-pg
LIBS = -lfftw3 -lm
HEADERS=Makefile headers/*.h

DIFFOPS=partials/partials.cpp findiff/findiff.cpp curl/curl.cpp div/div.cpp grad/grad.cpp lapl/lapl.cpp
UTILS=utility/bundlememcpy.cpp utility/linspace.cpp utility/pbc.cpp utility/printfield.cpp
MAINSRC=fftw_1d_test1.cpp
SOURCES=$(MAINSRC) $(DIFFOPS) $(UTILS)

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(HEADERS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o: $(HEADERS)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o *.exe
