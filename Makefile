CC = g++
C_FLAGS = -g

OBJECTS = main.o prototype.o

all: $(OBJECTS)
	$(CC) $(OBJECTS) -o main

main_icecube: main_icecube.o prototype.o
	$(CC) main_icecube.o prototype.o -o main_icecube

nonunitary: nonunitary.o prototype.o
	$(CC) nonunitary.o prototype.o -o nonunitary

%.o: %.cc %.h Makefile
	$(CC) $(C_FLAGS) -c $< -o $*.o

clean:
	-rm main nonunitary
	-rm *.o 

