CC=g++
CFLAGS=-c -Wall -std=c++14 -Iheaders -O3 #-pg
LDFLAGS= -O3 #-pg
LIBS = -lfftw3 -lm
HEADERS=Makefile headers/*.h

DIFFOPS=partials/partials.cpp findiff/findiff.cpp curl/curl.cpp div/div.cpp grad/grad.cpp lapl/lapl.cpp
UTILS=utility/fftwmemcpy.cpp utility/linspace.cpp utility/pbc.cpp utility/printfield.cpp utility/bundlememcpy.cpp
IO=io/writetofile.cpp io/savestepfname.cpp
SOLVERS=solvers/poisson_solver.cpp
MAINSRC=bundle_cpuderiv.cpp #solvers/kolmogorov_flow.cpp #solvers/kolmogorov_flow_debug.cpp
COMPARESRC=solvers/kflow_noPnclOptim_debug.cpp
SOURCES=$(MAINSRC) $(DIFFOPS) $(UTILS) $(IO) $(SOLVERS)
COMPARESOURCES=$(COMPARESRC) $(DIFFOPS) $(UTILS) $(IO) $(SOLVERS)

OBJECTS=$(SOURCES:.cpp=.o)
COMPOBJECTS=$(COMPARESOURCES:.cpp=.o)
EXECUTABLE=test.exe
COMPAREXEC=compare.exe

all: $(SOURCES) $(EXECUTABLE)

main: $(SOURCES) $(EXECUTABLE)

compare: $(COMPARESOURCES) $(COMPAREXEC)

$(EXECUTABLE): $(OBJECTS) $(HEADERS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@

$(COMPAREXEC): $(COMPOBJECTS) $(HEADERS)
	$(CC) $(LDFLAGS) $(COMPOBJECTS) $(LIBS) -o $@

.cpp.o: $(HEADERS)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o *.exe

cleanse:
	rm -f partials/*.o findiff/*.o curl/*.o div/*.o grad/*.o lapl/*.o utility/*.o io/*.o solvers/*.o


