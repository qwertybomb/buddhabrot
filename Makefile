CC = g++

FLAGS = -s -Ofast -march=native -fopenmp -std=gnu++20

LIBS = -lpng

EXE = main

buddhabrot: main.cc
	$(CC) $(FLAGS) main.cc -o $(EXE) $(LIBS)
clean:
	rm $(EXE)
