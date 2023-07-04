.. _output:

Output Files
============

Basics
------

In SCONE, the output facility for the calculation results is implemented through an *outputFile* class
that consists of the two main components. The first is the *outputClass* itself which checks the
sequence of output calls for consistency. The other is the printer (*asciiOutput_inter*), which
translates the output to a usable format. This division was done in order to enable writing output
into multiple formats through a common interface.

The structure of an output file is hierarchical. Its main component is the **block**, which is
similar to a dictionary. Each block is composed from a number of entries. Each is a pair of a name
together with some content. The content may be a **result**, **value**, **N-D array of results** or
**N-D array of values**. Value can be a **real**, **integer** or **character string**. Each result is
a pair of reals that represent a mean value of an estimate and its associated statistical uncertainty.
The initial block of output is called a *root block* and is not associated with any name in the output file.

.. admonition:: Note

  Relative vs Absolute standard deviation
    By convention, in SCONE, the uncertainty of statistical estimates in output files should be
    given as *absolute* values (e.g. ``4 +/- 2`` not ``4 +/- 50%``). As with any convention, it may be
    broken in exceptional circumstances, which must be *clearly* stated in the documentation of a
    relevant tallyClerk.

  Repeated Entry Names
    Names of each entry in a block must be *unique*. However, if a repeated name is present, only
    a warning is raised. Thus, an accidental name repetition will not kill the calculation with
    a fatal error (which would loose all results) and an output file will be produced. Although,
    it may parsed incorrectly by MATLAB, Python or other post-processing tool.

Writing to an output file in SCONE is done through a series of calls to appropriate procedures e.g.::

  character(nameLen) :: name
  type(outputFile)   :: out

  !! Initialise output by choosing the format and output file
  !! If no file name is provided, the output will not be printed
  call out % init('asciiMatlab', filename="./path/to/output/without_any_extension")

  !! Write an integer of 7
  name = 'int'
  call out % printValue(7, name)

  !! Long integers are fine
  name = 'large_int'
  call out % printValue(4294967296, name)

  !! We can open a nested block
  name = 'subBlock'
  call out % startBlock(name)

  !! And start an array
  name = 'myArray'
  call out % startArray(name, [3,2])

  !! The type of the array is determined by the first value provided
  !! Array is filled in column major order i.e.:
  !! 1st call -> index (1,1)
  !! 2nd call -> index(2,1)
  !! ...
  call out % addResult(7.8_defReal, 0.3_defReal)
  call out % addResult(7.8_defReal, 0.1_defReal)
  call out % addResult(7.8_defReal, 0.3_defReal)

  !! We can provide multiple numbers as an array for convenience
  call out % addResult([1.1_defReal, 1.2_defReal, 3.3_defReal], &
                       [0.1_defReal, 0.2_defReal, 0.3_defReal])

  !! We need to explicitly close array after all elements
  call out % endArray()

  !! Same goes for blocks
  call out % endBlock()

  !! We need to make sure that output file is properly finalised. We can skip
  !! this line if we we are not in a `program` block e.g. in `subroutine` or `function`
  call out % finalise()


The need to write name to a `character(nameLen)` variable before calling the procedure is a quirk
caused by the fact that the output file expects character of `nameLen` as a variable. If one were to
call the procedure with a shorter string literal (e.g. "subBlock") the missing characters would
be filled by a random memory content instead of blanks. This issue could be fixed by accepting the
deferred length character as an argument `character(*)`, which should be done in near future.

Error Checking
--------------

Note that the `outputFile` checks the correctness of a calls. So, for example, closing a block before
closing an array will produce a fatal error. The same will happen if type of entries in an array
changes between `add` calls. ::

  name = 'block'
  call out % startBlock(name)
  name = 'array'
  call out % startArray(name [2,3])

  !! Fatal error here
  call out % endBlock()

Similarly ::

  name = 'array'
  call out % startArray(name [5,3])
  call out % addValue(8.8_defReal)

  !! Fatal error here as well
  call out % addValue(7)

Other incorrect sequences of calls are also captured. Please refer to in-source documentation
of the `outputFile_class` for further details.

However, in many circumstances (e.g. unit tests) we may want to suppress the fatal error in case
of an incorrect output sequence. This can be done at initialisation by an optional
argument ::

  call out % init('asciiJSON', fatalErrors=.false.)

  ...
  ...
  Print output
  ...
  ...

  !! Check if there were errors
  if (.not.out % isValid()) then
    !! Print the error messages
    print *, out % getErrorLog()
  end if


MATLAB Format
-------------

Matlab format of SCONE produces a ``.m`` file with executable MATLAB code that defines a number of
variables. To achieve the hierarchical structure name prefixes are used. Thus, for example, an entry
'int' in block 'subBlock' and 'array' in root block will be printed as ::

   subBlock_int = 7;
   array = reshape([1,2,3,4,5,6],3,2);

The multidimensional arrays are printed as a 1D array and the final shape is achieved by the
`reshape` function. The results are printed as an array of two elements. Similarly an array
of results of dimension N will be printed as a N+1 dimensional array, of which the 1st column
(leftmost index ) will have length 2 and store the mean and the standard deviation.

Note that the output file is not intended to be read by a human in spite of being written in ASCII
characters. For example, every array is printed on a single line irrespective of its length, so
a user might find it rather difficult to inspect. The better way to read the output is to read it
into MATLAB by running it as a script ::

  run 'output.m'

  % Read standard deviation of a (1,4) entry in a result array
  res = block1_myResArray(2,1,4)



JSON Format
-----------

JSON format is intended as an arguably better alternative to the MATLAB output, that provides not
only better support for the hierarchy inherent to the output file, but also can be read without any
commercial third-party tool. In particular it can be easily imported into Python and processed with
NumPy. Output files that use JSON are produced with ``.json`` extension.

To read the JSON output file with Python the following code can be used ::

  import json
  import numpy as np

  with open('output.json') as f:
    # Reads output file into Python dictionary
    data = json.load(f)

  # Get an array from a block and print the standard deviation of (1,4) entry in a result array
  # Note that the transpose is required to recover the same indexing as in MATLAB
  Res = np.array(data['block1']['myResArray']).T

  # Avoid off-by-1 mistakes due to 0-indexing in Python
  print(Res[1, 0, 3])

The treatment of results and array of results in the output file is the same as in the MATLAB format.
However, when reading data into NumPy, please note that the order of indexes will be inverted,
as a result of the row-major order used in NumPy. Thus, it is necessary to take the transpose for
the indexing of a ``numpy.array`` to match the index order in MATLAB or supplementary information
provided by a tallyClerk.

When working with JSON output it is worth to consider using the `Spyder <https://www.spyder-ide.org/>`_
environment. In particular its variable explorer provides with an easy way to navigate through
the hierarchy of an output file.
