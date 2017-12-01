# Define compilation variables 
FC=gfortran
OBJ = ./../../objFiles
MOD = ./../../modFiles
BIN = ./../../binFiles 
FLAGS = -J $(MOD) -std=f2003 -c  

all: ByIsoNoMT_Data_class.o

ByIsoNoMT_Data_class.o:
	$(FC) $(FLAGS) ByIsoNoMT_Data_class.f03 -o $(OBJ)/ByIsoNoMT_Data_class.o
