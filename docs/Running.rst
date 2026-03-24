.. _running:

Running
=======

Note that instructions on how to write an input file can be found
in the :ref:`Input Manual<input-manual>`.

Running SCONE
'''''''''''''

After installation of all the dependencies and the compilation of SCONE, 
one can finally run the code. The executable ``scone.out`` is in the ``build``
folder; one can run it from any folder by specifying the path to the executables
and to the input file::

   <path_to_build_folder>/scone.out <input_file>

.. admonition:: Options

   OpenMP
     Specifies the number ``<n>`` of OpenMP threads to be used::

	--omp <n>

   Geometry plotting
     Allows plotting the geometry without running the actual calculations::

	--plot

   MPI
     To run with multiple processes, one needs to run using ``mpirun`` and 
     specify the number of processes <n>::

	<path_to_mpirun>/mpirun -np <n> <path_to_build_folder>/scone.out <input_file>

     OpenMP and MPI can be run together too::

	<path_to_mpirun>/mpirun -np <n> <path_to_build_folder>/scone.out --omp <n> <input_file>
