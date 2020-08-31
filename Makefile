CC = g++

OBJECTS = main.o prototype.o

all: $(OBJECTS)
	$(CC) $(OBJECTS) -o main

.c.o:
	$(CC) $(C_FLAGS) -c $< -o $*.o

clean:
	-rm main
	-rm *.o 

