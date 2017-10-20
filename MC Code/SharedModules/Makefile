# Define compilation variables 
FC=gfortran
OBJ = ./../objFiles
MOD = ./../modFiles
BIN = ./../binFiles 
FLAGS = -J $(MOD) -c  

all: numPrecision.o genericProcedures.o 

numPrecision.o:
	$(FC) $(FLAGS) numPrecision.f03 -o $(OBJ)/numPrecision.o

genericProcedures.o:
	$(FC) $(FLAGS) genericProcedures.f03 -o $(OBJ)/genericProcedures.o