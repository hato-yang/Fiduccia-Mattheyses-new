##############################################################
#Main compile for Fiduccia-Mattheyses
##############################################################
CFLAGS=-g
CC=gcc

all: src/main.o src/dll_structure.o src/basic_objects.o src/populate_partitions.o src/data_input.o src/fiduccia_mattheyses.o src/genetic_algorithm.o
	$(CC) $(CFLAGS) -o src/GAmain.out src/main.o src/dll_structure.o src/basic_objects.o src/populate_partitions.o src/data_input.o src/fiduccia_mattheyses.o src/genetic_algorithm.o
#	rm src/*.o


run: all
	src/GAmain.out

main: src/main.c include/main.h
	gcc  -c -g main.c include/main.h

dll_structure: src/dll_structure.c 
	gcc  -c -g dll_structure.c

basic_objects: src/basic_objects.c 
	gcc  -c -g basic_objects.c 

populate_partitions: src/populate_partitions.c
	gcc  -c -g populate_partitions.c populate_partitions.h

data_input: src/data_input.c 
	gcc  -c -g data_input.c 

fiduccia_mattheyses: src/fiduccia_mattheyses.c fiduccia_mattheyses.h
	gcc  -c -g fiduccia_mattheyses.c fiduccia_mattheyses.h

genetic_algorithm: src/genetic_algorithm.c genetic_algorithm.h
	gcc  -c -g genetic_algorithm.c genetic_algorithm.h

##############################################################
#Extra options
##############################################################
#Check heap memory for leaks
valgrind: all
	valgrind --leak -check=full --track-origins=yes src/main.out

