CXX = g++
CXXFLAGS = -O2 -march=native -std=c++17 -w
CC = gcc
CFLAGS = -O2 -w
SKETCHSRC = ArgusPS.cpp ArgusPH.cpp ArgusPP.cpp PSketch.cpp StableSketch.cpp Pandora.cpp WavingSketch.cpp Hypersistent.cpp OnOff.cpp PSSketch.cpp TightSketchPH.cpp
HEADERS = $(wildcard *.hpp) hash.h util.h

all: bench gen_synth

hash.o: hash.c hash.h
	$(CXX) $(CXXFLAGS) -c hash.c -o hash.o

bench: bench.cpp $(SKETCHSRC) $(HEADERS) hash.o
	$(CXX) $(CXXFLAGS) -o bench bench.cpp $(SKETCHSRC) hash.o

gen_synth: gen_synth.cpp
	$(CXX) $(CXXFLAGS) -o gen_synth gen_synth.cpp

clean:
	rm -f bench gen_synth hash.o
