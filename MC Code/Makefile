# Define compilation variables 
FC=gfortran
OBJ = ./objFiles
MOD = ./modFiles
BIN = ./bin
FLAGS = -J $(MOD) -c  

all: 
	make -C ./SharedModules
	make -C ./RandomNumbers
	make -C ./ParticleObjects
	make -C ./NuclearData

test: all
	$(FC) $(FLAGS) test.f03 -o $(OBJ)/test.o 
	$(FC) -J $(MOD) ./objFiles/*.o 
clean: 
	rm ./objFiles/*.o
	rm ./modFiles/*.mod
	rm ./*.out
	rm ./bin/*.out