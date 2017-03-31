BIN = ./bin/md-mpi
CC = mpicc
CFLAGS = -std=c99 -Wall -g 
INC = -I ./src 
SRC = $(wildcard src/*.c)
LIB = -lm

all: $(BIN)

$(BIN):$(SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(INC) $(LIB)

.PHONY:clean
clean:
	rm -rf $(BIN)
