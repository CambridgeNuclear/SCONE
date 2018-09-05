*****
SCONE
*****

SCONE (**S**\ tochastic **C**\ alculations **O**\ f **N**\ eutron Transport **E**\ quation) 
is an attempt to create an Object-Oriented framework for Monte Carlo particle transport 
calculations. Its primary goal is to support neutron transport problems in nuclear reactor 
physics. Object-oriented structure was selected to allow for better modularity and 
(hopefully) ability to modify solution algorithms with virtually no knowledge about 
the implementation details of large chunks of functionality (tallies, geometry, nuclear data 
handling). Fortran 2008 was chosen over other languages (mostly C++) for its performance combined 
with relative simplicity. 


Getting Started
===============
        
Prerequisites
-------------

* Modern Fortran Compiler 

  * gfortran 6.3 or higher 
* CMake 2.8.11.2 or higher 

------------
Installation
------------
1. Download the repositry. Run the following commands:: 

     git clone https://Mikolaj_Adam_Kowalski@bitbucket.org/Mikolaj_Adam_Kowalski/cued-mc-code.git  
    
2. Create build folder in the project directory(e.g. Debug):: 

     cd ./cued-mc-code
     mkdir Debug
   
3. Generate makefile with CMake and compile::

     cmake -E chdir ./Debug cmake ./..
     make -C Debug












