CC = g++

FLAGS = -s -Ofast -march=native -fopenmp -std=gnu++20

LIBS = -lpng

INCLUDE = -I./png++/

EXE = main

buddhabrot: main.cc
	$(CC) $(INCLUDE) $(FLAGS) main.cc -o $(EXE) $(LIBS)
clean:
	rm $(EXE)
