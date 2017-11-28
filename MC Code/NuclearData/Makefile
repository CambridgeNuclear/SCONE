# Define compilation variables 
FC=gfortran
OBJ = ./../objFiles
MOD = ./../modFiles
BIN = ./../binFiles 
FLAGS = -J $(MOD) -c  

all: nuclearData_class.o

nuclearData_class.o:
	$(FC) $(FLAGS) nuclearData_class.f03 -o $(OBJ)/nuclearData_class.o
