CC=icpc
CFLAGS=-c -O3 -Wall -std=c++11
LDFLAGS=-O3
SOURCES=src/Check.cc src/Chipmunk.cc src/Gauss_Legendre.cc src/Homogeneous_Problem.cc src/Parser.cc src/Sn_Stochastic.cc src/Sn_Transport.cc 
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=chipmunk

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm src/*.o
