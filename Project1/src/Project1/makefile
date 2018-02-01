CC=g++
CFLAGS= -std=c++14 -o3#-Wall -Wextra -Wpedantic -std=c11 -g
DIST = simulation.h system.h potential.h particle.h trialfunction.h Potentials/harmonicoscillator.h Systems/randomsystem.h DataDump/datadump.h Parameters/parameters.h
SOURCES=main.cpp simulation.cpp system.cpp particle.cpp potential.cpp trialfunction.cpp Potentials/harmonicoscillator.cpp Systems/randomsystem.cpp DataDump/datadump.cpp Parameters/parameters.cpp
BINS=simulation

all: simulation

run: simulation
	./simulation

simulation:
	$(CC) $(CFLAGS) $(SOURCES) $^ -o $@

clean:
	rm -f $(BINS)
