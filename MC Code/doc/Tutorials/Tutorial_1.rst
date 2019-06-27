
First Steps - Compiling new files and tests
===========================================

In this tutorial you will learn: 

* How to create new branch in Git 
* How to add new source file to compile with SCONE 
* How to create and compile a simple unit test for your code 
* How to add new files to Git and commit changes 

This tutorial assumes that you use a UNIX-like enviroment or Eclipse IDE. 
Moreover some knowlage of basic UNIX command line commands is essential.  

Creating new Git Branch
-----------------------

#. Let us first navigate to a main source code directory. This is the one that contains all folders like
   *Apps*, *DataStructures*, *SharedModules* etc. Type ``ls`` command to verify that you are in correct 
   Folder. 
   
#. We will proceed to make some changes. However, we also want to preserve an image of the code in its
   current, clean form. This is where version control with git comes into play. Into your console 
   type ``git branch`` and you should see somthing like this :: 
     
     dancoffBellTally
     errors-module
     makeDictGreatAgain
     * master
     missingENDFLaws
     newNuclearData        
 
#. With this command we have asked git to list all local branches (that is independent images of the code 
   repository, changes in one branch do not influance a code in another unless we specifically combine 
   some of them with ``git merge``). The ``*`` indicates that we are currently in a ``master`` branch. 
   We want to create a new branch for our tutorial. We can do it by creating a new branch with 
   ``git branch tutorials``, where ``tutorials`` could be replaced with any name for new branch. We 
   also need to switch to the new branch we just created. Do it with ``git checkout tutorials``. If you are 
   using Eclipse it is easier to manage branches with GUI. Right click on SCONE project in ``Project Explorer`` 
   and go to ``Team -> Switch to -> New branch``. Enter new branch name in the pop-up windows and 
   accept default settings. Now when you type ``git branch`` you should see somthing like this :: 
   
     dancoffBellTally
     errors-module
     makeDictGreatAgain
     master
     * tutorials
     missingENDFLaws
     newNuclearData               

Creating new Fortran module
---------------------------

#. Now that we have created space where we can easily seperate our modifications from the original 
   code, let us do some changes. First we will create a new directory for the files related to the 
   tutorial. Type ``mkdir TutorialFiles`` to create new folder. In Eclipse right-click on SCONE project 
   in ``Project Explorer`` and select ``New -> Folder``.  
   
#. Now we shall create a new source file. First navigate to the new Folder by ``cd TutorialFiles`` and 
   create new file with ``echo '' >> tutorial1_func.f90``. In Eclipse right-clik on the ``TutorialFiles`` 
   folder and select ``New-> Fortran Source File``. In pop-up window under *template* select ``<None>``. 
   One the file is created open it with your text editor of choise and fill it with content :: 
     
     module tutorial1_func 
       
       use numPrecision 
    
       implicit none 
       private 
       
       public :: hello 
       public :: quest
       
     contains 
     
       !!
       !! Prints nonsense text to the console 
       !!
       !! Args: 
       !!   None 
       !! 
       !! Errors: 
       !!   None 
       !!
       subroutine hello() 
        
         print *, "I wish for some good Devonshire cream tea..." 
       
       end subroutine hello 
     
       !!
       !! Return an integer 
       !!
       !! Args: 
       !!   None 
       !!
       !! Result: 
       !!   Answer to the Ultimate Question. 
       !!
       !! Errors: 
       !!   None  
       !!
       elemental function quest() result(i) 
         integer(shortInt)     :: i 
         
         i = 42 
         
       end function quest 
     
     end module tutorial1_func    

#. OK. So now we have a new Fortran module. Note that we are already following the Style Guide (:ref:`style-guide`). 
   In tutorial it is not necessary, but from personal experiance I recommend to try to get into 
   the habit of writing documentation and descriptive comments from the get-go. What we want to do now 
   is to compile the new module. To do that we need to register it with cmake so it is included in 
   a long list of source files that get compiled into libscone.a library. In the TutorialFiles folder 
   create new file ``CMakeLists.txt``. Note that the capitalisation and extension is important! 
   Inside the file write :: 
   
     # Add source files to a global list 
     add_sources(./tutorial1_func.f90)        

#. What we have done is we have created a file for CMake with instructions on what do to in this folder. 
   However, we haven't told CMake that ``TutorialFiles`` folder exists at all! To do this we need to step 
   to the folder below (``cd ..``), which is our case is the main source directory. Within it we can 
   easily find ``CMakeLists.txt`` file. In it we need to add a line ``add_subdirectory(TutorialFiles)``. 
   You can easily find the section where the other folders are registered :: 
   
      ###############################################################################
      # COLLECT ALL SOURCE AND TEST FILES

      # Include Nested Directories  
      add_subdirectory(RandomNumbers) 
      add_subdirectory(LinearAlgebra)
      add_subdirectory(SharedModules)    
      add_subdirectory(VTK)  
      add_subdirectory(ParticleObjects)  
      add_subdirectory(NamedGrids)

      add_subdirectory(NuclearData)      
      add_subdirectory(GeometryObjects)  
      add_subdirectory(Tallies)

      add_subdirectory(CollisionOperator)
      add_subdirectory(TransportOperator)

      add_subdirectory(MOCObjects) 
      add_subdirectory(UserInterface)

      add_subdirectory(PhysicsPackages)
      add_subdirectory(DataStructures)
      
      # Tutorial Folders  
      add_subdirectory(TutorialFiles)       
      ###############################################################################
      
#. Now your module should be compiled. Test that this is a case by introducing a deliberate error 
   in ``tutorial1_func.f90``.
   
   
Creating a simple Unit Test
---------------------------

#. Now that we have written some code we need to make sure that it does what we want. In practice 
   no matter how simple the function is and how unlikely it is that we have make a mistake some 
   sneaky bug will find its way in. It is best to catch it quickly! Then we will spend much less time 
   debugging some mysterious crashes. If you are not convinced by the need for automated testing yet 
   please have a look at extra arguments (:ref:`unit-testing`). 
   
#. First thing we need to do is to create folder to store our test files. In ``TutorialFiles`` create 
   ``Tests`` subfolder. Inside it create a ``tutorial1_test.f90`` file. And fill it as follows :: 
   
     module tutorial1_test
      
       use tutorial1_func, only : quest  
       use pFUnit_mod 
     
       implicit none 
     
     contains 
       
       !!
       !! Test the Ultimate Question 
       !!
     @Test 
       subroutine testUltimateQuestion() 
       
         !! Test Question 
         @assertEqual( 43, quest() )
         
       
       end subroutine testUltimateQuestion 
     
     
     end module tutorial1_test   

#. You may notice that the answer to the question is 42 not 43, so the test will fail for the correct 
   result. This is deliberate. We want the test to fail so we can make sure that the test is indeed 
   compiling and executing. Then we can change the values back to what we expect and ensure that the 
   code is indeed correct. In order to finish we must register the test with out CMake script. We need
   to go back to the ``CMakeLists.txt`` in the ``TutorialFiles`` folder and add a line :: 
   
     # Add source files to a global list 
     add_sources(./tutorial1_func.f90)   
     
     # Register unit tests 
     add_unit_tests(./Tests/tutorial1_test.f90) 
     
#. Now recompile the code making sure that you are compiling it together with tests! Upon sucessfull 
   compilation execute ``allTests`` binary in the build folder. It should give you a message like 
   this :: 
   
     ................................................................................
     ...........................................................F
     Time:         0.079 seconds
  
     Failure
      in: 
     tutorial1_test_suite.testUltimateQuestion
      Location: 
     [tutorial1_test.f90:11]
     expected 43 but found: 42;  difference: |1|.
  
      FAILURES!!!
     Tests run: 139, Failures: 1, Errors: 0

                
#. Clearly our new test was run and failed! Now go back to ``tutorial1_test.f90`` and change 43 to 42. 
   Recompile and execute the test binary again. Now it should read :: 
   
     ................................................................................
     ...........................................................
     Time:         0.101 seconds
  
     OK
     (139 tests)

#. Which means that the test was sucesfull and our trivial procedure does what we have expected. 
   Note that we have not tested ``hello`` subroutine. It is becouse it is immpossible to test the 
   console messages within the pFUnit test framework. The same, unfortunately, goes for any 
   warnings or fatal errors that can be produced. Go here **ADD LINK** for more details. 
   
Creating an executable 
---------------------- 

#. Often we want to be able to execute some simple code to test that its behaviour is as we belive it 
   to be, or we want to play a bit with some of the SCONE components to get used to them. Using the 
   pFUnit test framework for it would be inconvinient. It is much simpler to create a simple program. 
   In order to do so we need to create a new file in Apps folder e.g. ``sandbox.f90``. :: 
     
     program sandbox 
     
       use numPrecision
       use tutorial1_func, only : hello  
       implicit none 
       
       call hello() 
       
     end program sandbox 
   
#. It is a bit more difficoult to register it with CMake and make sure it is compiled. Once again 
   go to the ``CMakeLists.txt`` in main source directory and add an extra entry to *COMPILE SOLVERS* 
   section :: 
   
     ###############################################################################
     # COMPILE SOLVERS 
     add_executable(test.out ./Apps/test_MAK.f90)
     target_link_libraries(test.out scone )  
     
     # TUTORIAL EXECUTABLE
     add_executable(sandbox.out ./Apps/sandbox.f90) 
     target_link_libraries(sandbox.out scone) 
     
     #add_executable(eigenCE.out        ./Apps/eigenCE.f90           ${SRCS}) 
     #add_executable(eigenMG.out        ./Apps/eigenMG.f90           ${SRCS}) 
     #add_executable(test_init_geom.out ./Apps/test_init_geom.f90    ${SRCS}) 
     #add_executable(dictTest.out       ./Apps/dictTest.f90          ${SRCS})
     #add_executable(parser.out         ./Apps/testACE.f90           ${SRCS})
     #add_executable(eigenvalue.out     ./Apps/eigenvalue.f90        ${SRCS})
     #add_executable(surfTests.f90       ./Apps/surfTests.f90      ${SRCS})
     #add_executable(geomPrint.out       ./Apps/geomPrint.f90      ${SRCS})
     #add_executable(rand.out        ./Apps/produceRNG.f90 ${SRCS})

     add_executable(scone.out ./Apps/scone.f90 )
     target_link_libraries(scone.out scone )       

     ###############################################################################   
   
#. Now recompile and a new executable ``sandbox.out`` should be produced in build folder. Run it 
   and see that it indeed produces the expected result.      
   
Commiting changes
-----------------
#. In the last step we would like to add the changes we have made to the git repository. We have 
   created a number of new files so we need to make Git track their changes. First of all let us 
   see what new files Git sees. Write ``git status`` and you should see something like this :: 
   
     # On branch tutorials
     # Changes not staged for commit:
     #   (use "git add <file>..." to update what will be committed)
     #   (use "git checkout -- <file>..." to discard changes in working directory)
     #
     #       modified:   CMakeLists.txt
     #
     # Untracked files:
     #   (use "git add <file>..." to include in what will be committed)
     #
     #       Apps/sandbox.f90
     #       Debug/
     #       TutorialFiles/
     no changes added to commit (use "git add" and/or "git commit -a")

#. So clearly Git sees that we have created new files. And modified some it is already tracking. You
   might have already heard of the term *commit*. Basicaly each commit is a snapshot of a state 
   of repository. Each branch is then a chain of commits (snapshots) from the latest one to the 
   first one. Thus to save our changes we need to create a new commit. We begin by *staging* all 
   the changes we have made. You can do it file-by-file using ``git add <file>`` as the git status 
   message suggests. Note that we cannot use ``git add --all`` becouse we have to ensure that Debug 
   folder that contains all compilation files in not tracked to keep repository clean.  Now when we 
   use ``git status`` again we shall see somthing like this:: 
       
     # On branch tutorials
     # Changes to be committed:
     #   (use "git reset HEAD <file>..." to unstage)
     #
     #       new file:   Apps/sandbox.f90
     #       modified:   CMakeLists.txt
     #       new file:   TutorialFiles/myFile
     #

#. Now what is left do do is commit the changes. We use ``git commit`` for that. If we type it 
   just like it is it will open a text editor in the console. Usually it will be vi. We need to 
   prepare commit message. Use *Insert* to enter edit mode in vi. In the first line provide a 
   basic explenation for the changes e.g. *Created Tutorial Infrastructure*. Make sure it is 
   shorter then 50 columns. If you want to provide same extra information, leave 2nd line blank 
   and write some extra text starting from the 3rd line. After you are finished use *Esc* key to 
   exit from vi edit mode and type ``:wq`` to write your changes and quit. After that you will 
   see somthing like that ::           

     [tutorials 5b8b5bd] Created some random files
     3 files changed, 4 insertions(+), 2 deletions(-)
     create mode 100644 MC Code/Apps/sandbox.f90
     create mode 100644 MC Code/TutorialFiles/myFile
  
#. You have created the commit. Note that if you want to provide only the basic message with your 
   commit you can do it without using vi with ``git commit -m '<Your message of less the 50 characters>'``. 
      
#. The entire process looks a lot diffrent with Eclipse IDE. First we need to add files one-by one 
   by right-cliking on them and going to ``Team -> Add to index``. After we are finished we can 
   create commit by right-cliking on SCONE project in Project Explorer and going to 
   ``Team -> Commit...``. Then in the pop-up windows we can select all the changes we want to include. 
   Once we are done we can just click ``commit``.
   
#. This concludes the tutorial.    