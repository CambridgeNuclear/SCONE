
Level 0 - Getting Started with Fortran 
======================================

In this tutorial you will learn:

* How to compile a simple Fortran program with gfortran compiler 
* Learn number of Fortran 2008 features used extensively in SCONE 

This tutorial assumes that you use a UNIX-like enviroment.
Moreover some knowlage of basic UNIX command line commands is essential. 
Furthermore gfortran in version 6.3 or higher should be installed on your system. 

Fortran 2008 Hello World 
------------------------

#. Let us start by crating a new file with ``echo '' >> helloFortran.f90``. Open the new file 
   in your editor of choice and write the following hello world program. :: 
     
     program helloFortran
       implicit none 
       
       print *, "Hello. I am Fortran. The first of high-level programing languges."  
     
     end program helloFortran 
  
#. Now verify the version of your compiler by typing ``gfortran --version``. You should see 
   somthing like this :: 
   
     GNU Fortran (GCC) 6.3.0
     Copyright (C) 2016 Free Software Foundation, Inc.
     This is free software; see the source for copying conditions.  There is NO
     warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.   
     
#. Now let us compile with following command ``gfortran helloFortran.f90 -o hello.out`` to create 
   new executable *hello.out*. Run it to verify that it produced the expected result. 
   
#. Now you should have no problem to learn Fortran basic syntax. For instance use this 
   `tutorial <https://www.tutorialspoint.com/fortran/>`_. 
   
Automatic Allocation 
--------------------     