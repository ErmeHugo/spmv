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
	#icc -o spmv spmv_mkl.c -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
	./$(PROG)
	$(RM) $(OBJ) $(PROG)

.PHONY : clean save

clean :
	$(RM) $(OBJ) $(PROG)
	
save : 
	mkdir save/
	cp -f $(SRC) save/
