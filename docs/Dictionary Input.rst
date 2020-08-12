.. _dictSyntax:

Dictionary Input
================

ASCII Dictionary Syntax
-----------------------

SCONE uses a hierarchical input that is composed of nested dictionaries
that contain some content associated with a keyword unique within the scope
of a single dictionary. Following content is available:

* Word with no spaces or ";" e.g. ThisIsContent
* Integer number e.g. 1
* Real number e.g. 2.78
* List of Words, Integers or Reals
* Subdictionary

Hierarchical structure of dictionaries can be loaded from ASCII files and be
used inside the SCONE to construct different objects and set calculation
settings. The syntax for writing dictionaries is based on
`OpenFOAM <https://cfd.direct/openfoam/user-guide/basic-file-format/>`_. However only subset of
OpenFOAM syntax is supported. An example of the correct input dictionary is::

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
      greeks{ poet Homer; politician Pericles; hero Theseus; }
      ! Note lack of semicolon at the end of subdictionary

At least one space is needed before an entry name and its value ::

    integer 1    ;    - OK
    integer        1; - OK
    integer1;         - WRONG
    intArray ( 1 2 3) - OK
    intArray(1 2 3)   - WRONG

Note that this is not a case for sub-dictionaries ::

    sub{ entry1 1;}   - OK
    sub { entry1 1;}  - OK

It is important to note that:

#. All TABS and NEWLINE are converted to single SPACE
#. NEWLINE should be UNIX. Convert file with dos2unix utility
#. File cannot contain any Unicode characters. e.g.: é Æ Λ
#. Input is free form but anything beyond column 1000 will be ignored!
#. Currently on 32bit integers can be read (max 2'147'483'647)
#. 64-bit Integer will be read as a character.
#. Maximum length of character in list is 30. Character not in list can have length of 100.
#. Real entry like ``70E+03`` will be read as real ``~0.0`` . Each real number must contain a dot
   ``70.E+03`` will be fine.

Some of these rules may change in near-future. Any suggestion for improvements are welcome.

Dictionary Grammar
------------------

It is not particularly useful by might be a good starting point for somebody sometime. A sketch
of SCONE dictionary grammar in BNF notation is::

    <dictionary>  ::= " "* <item>+ " "*
    <item>        ::= <entry> | <word> " "* "{" <dictionary> "}"
    <entry>       ::= <word> " "+ <content> " "* ";"
    <content>     ::= <single> | "(" <list> ")"
    <single>      ::= <int> | <real> | <word>
    <list>        ::= <intList> | <realList> | <wordList>
    <intList>     ::= (" "* <int> " "+)+
    <realList>    ::= (" "* (<int>|<real>)" "+)* (" "* <real>" "+)+ (" "* (<int>|<real>)" "+)*
    <wordList>    ::= (" "* <word> " "+)+

For compactness definitions of ``<int>``, ``<real>`` and ``<word>`` are omitted.
They have the following meaning::

    <int> ::= Characters representing and integer e.g.: ^\d+$
    <real> ::= Real number e.g.: ^[-+]?[0-9]+[.][0-9]*([eE][-+]?[0-9]+)?$
    <word> ::= Must contain not digit character e.g.: .*[a-zA-Z]+.*

This is not yet full definition of the grammar as it does not contain limitations related to
maximum size of integers, reals and maximum length of the characters. It may be useful if somebody
will try to use a proper parser to read the SCONE dictionary files.
