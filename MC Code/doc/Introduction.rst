
Introduction
============

What is SCONE
-------------

SCONE stands for **S**\ tochastic **C**\ alculator **O**\ f **N**\ eutron Transport **E**\ quation and  
it is a Monte Carlo particle transport code, like its more famous cousins Serpent, 
MCNP, OpenMC ,MCS, MONK, Tripoli, MORET, RMC and many others. Naturally, the existance of so many good and 
established Monte Carlo codes can make anybody question the need for **yet** another code. 
However, most of the MC codes are proprietary and access to them is quite limited. Moreover most 
of them strive for maximum performance, which can reduce a readibilty of their source code for 
a novice user unexperianced with programming. In our view, this creates a niche for a highly flexible 
and easy-to-modify MC neutron transport code, that would allow very novice users like Master and 
PhD students to test new approches to preforming Monte Carlo calculations, even if they have no 
privuous significant programming background. This is of paramount importance in the context of British 
Masters programmes that last only a single year, so it is necessary a for student to get a grasp of 
their programming enviroment quickly in order to be able to do any significant work on their project. 

With this in mind aims of SCONE are as follows: 

* Minimise time between an idea and prototype implementation on semi-realistic test case. 
* Be easy to learn from programming perspective.  
* Provide detailed explenation of algorithims and physics.  
* Have easy to modify source code, even if it makes input files more complicated. 
* Have reasonable, not maximum performance. Trade-offs between modifiability and performance should 
  be skewed towards modyfiability.  

Overall, it is perhaps better to consider SCONE an advanced "toy code" rather than a software for 
real calculations. The day SCONE will be properly validated will be the day a brand new AGR 
will be build!   

For the british people out there, we would like to stress that we do not have a prefered pronunciation
both the *tone* and *gone* rhyme is perfectly fine. It is however a strong opinion of this author that 
the jam goes on top of clotted cream!     

Installation Instructions
-------------------------

Requirements
''''''''''''

Following is required to compile SCONE: 

* Cmake version 3.0 or higher 
* Modern Fortran Compiler (gfortran 6.3 or higher)
* LAPACK and BLAS Libraries 
* pFUnit unit testing framework and Python interpreter (at least 2.7)  
* UNIX-like enviroment. Linux, MacOS or Cygwin in Windows 

Getting Fortran Compiler
''''''''''''''''''''''''
If you happen to work on a computer that does not have an up-to-date Fortran compiler. You need to 
obtain it. 

* **Add Instalation instructions with root access** 

However, if you don't have superuser privileges you may still compile the GCC from source: 

* **Add Instructions here**


Installing pFUnit
'''''''''''''''''
This is only required if the unit tests are to be build. 

#. Make sure python can be invoked by a command 'python' by typing:: 

     python --version 

#. Create a folder for the local installation of pFUnit e.g. in your home directory and 
   download the pFUnit repository and enter the source code folder:: 
   
     mkdir pFUnit
     cd pFUnit
     git clone git://git.code.sf.net/p/pfunit/code pfunit-code
     cd pfunit-code
          
#. Export environmental variables required by pFUnit:: 

     export F90=gfortran
     export F90_VENDOR=GNU  
     
#. Build and test pFUnit by typing::

     make tests 
     
#. Install pFUnit in any directory you have access to e.g. :: 

     make install INSTALL_DIR=~/pFUnit
     
LAPACK and BLAS
'''''''''''''''
#. Download a version of LAPACK from `official website <http://www.netlib.org/lapack/>`_.

#. In some directory on your filesystem extract the archive.

#. Configure compilation with cmake by typing:: 

     mkdir Build 
     cd Build
     cmake ./..

#. If you don't have a root access on your machine or you want to install LAPACK to  a custom 
   directory, use ccmake to change CMAKE_INSTALL_PREFIX. In Build directory type::
   
     ccmake ./..  
     <Navigate to CMAKE_INSTALL_PREFIX and change it to your folder> 
     Press [c] to configure 
     Press [g] to generate and exit 
     
#. Now compile LAPACK and install by typing:: 

     make 
     make install      
     
     
Compiling SCONE
'''''''''''''''

#. If you want to install with tests set PFUNIT_INSTALL environmental variable to directory in 
   which pFUnit was installed. e.g. :: 
   
     export PFUNIT_INSTALL=~/pFUnit    

#. If your LAPACK installation is not in default system directories use LAPACK_INSTALL enviromental 
   variable to help CMAKE find the library. e.g. :: 
   
     export LAPACK_INSTALL=~/LAPACK 

#. Download the repositry. Run the following commands:: 

     git clone https://Mikolaj_Adam_Kowalski@bitbucket.org/Mikolaj_Adam_Kowalski/scone.git  
    
#. Create build folder in the project directory(e.g. Debug):: 

     cd ./cued-mc-code
     mkdir Debug
   
#. Generate makefile with CMake and compile::

     cmake -E chdir ./Debug cmake ./..
     make -C Debug

#. To switch off compilation of tests use the following commands:: 

     cmake -E chdir ./Debug cmake ./.. -DBUILD_TESTS=OFF 
     make -C Debug 

#. Note that you can use ccmake utility to modify avalible options and regenerate your make file just 
   type the following into your terminal and follow the instructions:: 

     ccmake ./Debug     

     
Getting Started
---------------
Having an Integrated Development Enviroment (IDE) can significantly help with the development. 
To develop SCONE we are using `Eclipse PTP <https://www.eclipse.org/ptp/downloads.php>`_. It contains 
a decent Fortran editor and GUI support for git (Egit). Unfortunatly Eclipse lack an integrated 
support for CMake, thus it is neccessary to create few workaround to use it.   

*  **Ad Instructions on Setting up a new project**