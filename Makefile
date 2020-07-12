CC = g++

FLAGS = -Ofast -march=native -mtune=native -fopenmp -std=c++20

LIBS += -lpng

INCLUDE += -I./png++/

EXE = main

buddhabrot: main.cc
	$(CC) $(INCLUDE) $(FLAGS) main.cc -o $(EXE) $(LIBS)
clean:
	rm $(EXE)