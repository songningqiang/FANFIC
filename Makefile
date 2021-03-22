CC = g++-10
C_FLAGS = -g -I/usr/local/Cellar/boost/1.74.0/include -std=c++11 -fopenmp

OBJECTS = main.o prototype.o

all: $(OBJECTS)
	$(CC) $(C_FLAGS) $(OBJECTS) -o main

nonunitarity: nonunitarity.o prototype.o
	$(CC) nonunitarity.o prototype.o -o nonunitarity

infersourcecomposition: infersourcecomposition.o prototype.o
	$(CC) infersourcecomposition.o prototype.o -o infersourcecomposition

infersourcefraction: infersourcefraction.o prototype.o
	$(CC) infersourcefraction.o prototype.o -o infersourcefraction

infersourcefraction_kpikmu: infersourcefraction_kpikmu.o prototype.o
	$(CC) infersourcefraction_kpikmu.o prototype.o -o infersourcefraction_kpikmu

neutrinodecay: neutrinodecay.o prototype.o
	$(CC) $(C_FLAGS) neutrinodecay.o prototype.o -o neutrinodecay

neutrinodecay_masseigenstates: neutrinodecay_masseigenstates.o prototype.o
	$(CC) $(C_FLAGS) neutrinodecay_masseigenstates.o prototype.o -o neutrinodecay_masseigenstates

neutrinodecay_kpigaussian: neutrinodecay_kpigaussian.o prototype.o
	$(CC) $(C_FLAGS) neutrinodecay_kpigaussian.o prototype.o -o neutrinodecay_kpigaussian


%.o: %.cc prototype.h Makefile
	$(CC) $(C_FLAGS) -c $< -o $*.o

clean:
	-rm *.o 
	-rm main nonunitarity infersourcecomposition infersourcefraction infersourcefraction_kpikmu \
	neutrinodecay neutrinodecay_masseigenstates  neutrinodecay_kpigaussian
	

