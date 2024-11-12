all : spmv

CC = gcc -Wall
RM = rm -f
SRC = $(wildcard *.c)
HEAD = $(wildcard *.h)
OBJ = $(SRC:.c=.o)
OPT = -O0
PROG= spmv


spmv : $(OBJ)
	$(CC) $(OPT) $^ -lopenblas -lgsl -lgslcblas -o $(PROG)
	./$(PROG)
	$(RM) $(OBJ) $(PROG)

.PHONY : clean save

clean :
	$(RM) $(OBJ) $(PROG)
	
save : 
	mkdir save/
	cp -f $(SRC) save/
