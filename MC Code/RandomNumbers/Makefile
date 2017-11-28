# Define compilation variables 
FC=gfortran
OBJ = ./../objFiles
MOD = ./../modFiles
BIN = ./../binFiles 
FLAGS = -J $(MOD) -c  

all: RNG_class.o

RNG_class.o:
	$(FC) $(FLAGS) RNG_class.f03 -o $(OBJ)/RNG_class.o
