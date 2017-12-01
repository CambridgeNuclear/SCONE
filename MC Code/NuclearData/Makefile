# Define compilation variables 
FC=gfortran
OBJ = ./../objFiles
MOD = ./../modFiles
BIN = ./../binFiles 
FLAGS = -J $(MOD) -c  

all: ByIsoNoMT

nuclearData_class.o:
	$(FC) $(FLAGS) nuclearData_class.f03 -o $(OBJ)/nuclearData_class.o

ByIsoNoMT: nuclearData_class.o
	make -C ./ByIsoNoMT