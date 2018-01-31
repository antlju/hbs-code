CC=g++
CFLAGS=-Wall -std=c++11 -O0
LDFLAGS=
DER1SOURCES=der1test_wdivergence_v2.cpp derivatives.cpp #Dependencies .cpp
ABCSRC=abctest_v2.cpp derivatives.cpp
CURLSRC=curldev_v1.cpp derivatives.cpp
HDEPS=common.h field.h
OBJECTS=$(SOURCES:.cpp=.o) #Dependencies .o
EXECUTABLE=der1test_wdivergence_v2.cpp

sinderivtest.o: derivsintest.cpp derivatives.cpp common.h field.h
	$(CC) $(CFLAGS) derivsintest.cpp derivatives.cpp -o $@

der1_divtest.o: $(DER1SOURCES) $(HDEPS)
	$(CC) $(CFLAGS) $(DER1SOURCES) -o $@

abctest_v2.o: $(ABCSRC) $(HDEPS)
	$(CC) $(CFLAGS) $(ABCSRC) -o $@

curldev_v1.o: $(CURLSRC) $(HDEPS)
	$(CC) $(CFLAGS) $(CURLSRC) -o $@

clean:
	rm *.o

