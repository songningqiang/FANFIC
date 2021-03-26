CC = g++-10
C_FLAGS = -g -I/usr/local/Cellar/boost/1.74.0/include -Iinclude -std=c++11 -fopenmp

SRC_DIR := src
OBJ_DIR := src
BIN_DIR := .

OBJECTS = $(OBJ_DIR)/prototype.o

EXES = main nonunitarity infersourcecomposition infersourcefraction infersourcefraction_kpikmu \
	neutrinodecay neutrinodecay_masseigenstates neutrinodecay_kpigaussian

all: $(EXES)

main : $(OBJ_DIR)/main.o $(OBJECTS) | $(BIN_DIR)
	$(CC) $(C_FLAGS) $^ -o $@

nonunitarity : $(OBJ_DIR)/nonunitarity.o $(OBJECTS) | $(BIN_DIR)
	$(CC) $(C_FLAGS) $^ -o $@

infersourcecomposition : $(OBJ_DIR)/infersourcecomposition.o $(OBJECTS) | $(BIN_DIR)
	$(CC) $(C_FLAGS) $^ -o $@

infersourcefraction : $(OBJ_DIR)/infersourcefraction.o $(OBJECTS) | $(BIN_DIR)
	$(CC) $(C_FLAGS) $^ -o $@

infersourcefraction_kpikmu : $(OBJ_DIR)/infersourcefraction_kpikmu.o $(OBJECTS) | $(BIN_DIR)
	$(CC) $(C_FLAGS) $^ -o $@

neutrinodecay : $(OBJ_DIR)/neutrinodecay.o $(OBJECTS) | $(BIN_DIR)
	$(CC) $(C_FLAGS) $^ -o $@

neutrinodecay_masseigenstates : $(OBJ_DIR)/neutrinodecay_masseigenstates.o $(OBJECTS) | $(BIN_DIR)
	$(CC) $(C_FLAGS) $^ -o $@

neutrinodecay_kpigaussian : $(OBJ_DIR)/neutrinodecay_kpigaussian.o $(OBJECTS) | $(BIN_DIR)
	$(CC) $(C_FLAGS) $^ -o $@


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc include/prototype.h Makefile
	$(CC) $(C_FLAGS) -c $< -o $@


clean:
	-rm $(OBJ_DIR)/*.o 
	-rm $(EXES)
	

