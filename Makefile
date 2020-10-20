EXEC := adjoint
SRC  := $(wildcard *.f90)
OBJ  := $(patsubst %.f90, %.o, $(SRC))
# NOTE - OBJ will not have the object files of c codes in it, this needs to be improved upon.
# Options	
F90 	:= gfortran
CC 		:= gcc
POP_PUSH:= ./pop_push

# Rules

$(EXEC): $(OBJ) adBuffer.o adStack.o
		$(F90) -o $@ $^

%.o: %.f90 
		$(F90) -c $<

driver.o: forward_b.f90 forward_diff.mod
forward_diff.mod: forward_b.o

adBuffer.o: 
		$(CC) -c $(POP_PUSH)/adBuffer.c
adStack.o : 
		$(CC) -c $(POP_PUSH)/adStack.c

forward_b.f90: forward.f90
		tapenade -reverse -head "forward_problem(V)/(xx)" forward.f90

# Useful phony targets

.PHONY: clean

clean:
	$(RM) $(EXEC) *.o $(MOD) $(MSG) *.msg *.mod *_b.f90