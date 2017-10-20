# Define compilation variables 
FC=gfortran
OBJ = ./../objFiles
MOD = ./../modFiles
BIN = ./../binFiles 
FLAGS = -J $(MOD) -c  

all: particle.o

particle.o:
	$(FC) $(FLAGS) particle.f03 -o $(OBJ)/particle.o
