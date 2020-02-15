.. _dictSyntax:

Dictionary Input
================

SCONE uses a hierarchical input that is composed of nested dictionaries
that contain some content associated with a keyword unique within the scope
of a single dictionary. Following content is available:

* Word with no spaces or ";" e.g. ThisIsContent
* Integer number e.g. 1
* Real number e.g. 2.78
* 1D Array of Words, Integers or Reals
* Subdictionary

Hierarchical structure of dictionaries can be loaded from ASCII files and be
used inside the SCONE to construct different objects and set calculation
settings. The syntax for writing dictionaries is based on OpenFOAM.::

      // This is a line comment. C++ Style
      ! Fortran Style line comments can also be used

      integer 1;   ! Semicolon separates entries
      intArr (1 2 3);
      real 4.3;
      realArr (1 2 3.4); ! Integers will be converted to reals

      int1 2; int2 3; ! Multiple entries can be in single line

      word Horace;
      words (Non omnis moriar);

      ! Subdictionary
      greeks { poet Homer; politician Pericles; hero Theseus; }
      ! Note lack of semicolon at the end of subdictionary

At least one space is needed before an entry name and its value ::

    integer 1    ;    - OK
    integer        1; - OK
    integer1;         - WRONG

During reading all tabs and newline characters are converted into a single space
Behaviour with DOS newline characters was not tested. Will probably fail.
Convert to unix with dos2unix utility!
