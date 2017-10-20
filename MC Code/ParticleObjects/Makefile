# Define compilation variables 
FC=gfortran
OBJ = ./../objFiles
MOD = ./../modFiles
BIN = ./../binFiles 
FLAGS = -J $(MOD) -c  

all: particle_class.o

particle_class.o:
	$(FC) $(FLAGS) particle_class.f03 -o $(OBJ)/particle_class.o
