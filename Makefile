all : spmv

CC = gcc -Wall
RM = rm -f
SRC = $(wildcard *.c)
HEAD = $(wildcard *.h)
OBJ = $(SRC:.c=.o)
PROG= spmv


spmv : $(OBJ)
	$(CC) $^ -lopenblas -lgsl -lgslcblas -o $(PROG)
	./$(PROG)
	$(RM) $(OBJ) $(PROG)

.PHONY : clean save

clean :
	$(RM) $(OBJ) $(PROG)
	
save : 
	mkdir save/
	cp -f $(SRC) save/
