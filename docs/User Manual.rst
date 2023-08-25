.. _user-manual:

User Manual
===========

Generic information about how to use dictionaries in writing an input file can be found
in :ref:`Dictionary Input <dictSyntax>`. Here, more specific information about the input
options available are described.

Physics Package
---------------

In SCONE, the **physics package** is what determines the flow of a calculation. The possible
options and related input files are shown below.

eigenPhysicsPackage
###################

eigenPhysicsPackage, used for criticality (or eigenvalue) calculations

* pop: number of particles used per cycle
* active: number of active cycles
* inactive: number of inactive cycles
* dataType: determines type of nuclear data used; can be ``ce`` or ``mg``
* XSdata: keyword to the name of the nuclearDataHandle used
* seed (*optional*): initial seed for the pseudo random number generator
* outputFile (*optional*, default = 'output'): name of the output file
* outputFormat (*optional*, default = ``asciiMATLAB``): type of output file.
  Choices are ``asciiMATLAB`` and ``asciiJSON``

Example: ::

        type eigenPhysicsPackage;
        pop    100000;
        active 100;
        inactive 50;
        dataType ce;
        XSdata   ceData;
        seed     -244654;
        outputFile PWR_1;
        outputFormat asciiJSON;

        transportOperator { <Transport operator definition> }
        collisionOperator { <Collision operator definition> }
        inactiveTally { <Inactive tally definition> }
        activeTally { <Active tally definition> }
        geometry { <Geometry definition> }
        nuclearData { <Nuclear data definition> }

*Optional entries* ::

        uniformFissionSites { <Uniform fission sites definition> }
        varianceReduction { <Weight windows definition> }
        source { <Source definition> }

.. note::
   Although a ``source`` definition is not required, it can be included to replace
   the default uniform fission source guess used in the first cycle

fixedSourcePhysicsPackage
#########################

fixedSourcePhysicsPackage, used for fixed source calculations

* pop: number of particles used per batch
* cycles: number of batches
* dataType: determines type of nuclear data used. Can be ``ce`` or ``mg``
* XSdata: keyword to the name of the nuclearDataHandle used
* seed (*optional*): initial seed for the pseudo random number generator
* outputFile (*optional*, default = 'output'): name of the output file
* outputFormat (*optional*, default = ``asciiMATLAB``): type of output file.
  Choices are ``asciiMATLAB`` and ``asciiJSON``

Example: ::

        type fixedSourcePhysicsPackage;
        pop    100000;
        cycles 200;
        dataType ce;
        XSdata   ceData;
        seed     2829741;
        outputFile shield_type11;

        transportOperator { <Transport operator definition> }
        collisionOperator { <Collision operator definition> }
        tally { <Tally definition> }
        source { <Source definition> }
        geometry { <Geometry definition> }
        nuclearData { <Nuclear data definition> }

*Optional entries* ::

        varianceReduction { <Weight windows definition> }

rayVolPhysicsPackage
####################

rayVolPhysicsPackage, used to perform ray-tracing based volume calculation

* pop: number of rays used per cycle
* cycles: number of cycles
* mfp: mean length of ray segments
* abs_prob: ray absorption probability after each segment
* robust: 1 for true; 0 for false; enable robust mode: in this case at each collision,
  each particle verifies that the material it currently thinks it is in and the one
  obtained by *placing* a particle in the geometry with the same spatial position and
  direction are in agreement
* cache: 1 for true; 0 for false; enable distance caching
* seed (*optional*): initial seed for the pseudo random number generator

Example: ::

        type rayVolPhysicsPackage;
        pop    1000000;
        cycles 100;
        mfp    0.3;
        abs_prob 0.1;
        robust   1;
        cache    1;

        geometry { <Geometry definition> }
        nuclearData { <Nuclear data definition. Requires material names only> }

vizPhysicsPackage
#################

vizPhysicsPackage, used for visualising geometry

Example: ::

        type vizPhysicsPackage;

        geometry { <Geometry definition> }
        viz { <Visualiser definition> }

Source
------

For the moment, the only possible external **source** type in SCONE in a point source.
The properties of a point source are:

* r: (x y z) vector with the origin position. [cm]
* particle: ``neutron`` or ``photon``, according to the type of particles emitted by the
  source
* E or G: emission energy

  - E: energy of the particles emitted, for continuous energy calculations. [MeV]
  - G: energy group of the particles emitted, for multi-group calculations

* dir (*optional*, default = isotropic): (u v w) vector with the direction of the source
  particles

Hence, an input would look like: ::

      source { type pointSource; r (0.0 1.0 5.2); particle neutron; E 14.1; dir (0.0 1.0 0.0); }

Transport Operator
------------------

The **transport operator** takes care of moving the particles from one collision location
to another. In the input file, one must include: ::

      transportOperator { type <transportOperatorType>; *keywords* }

The possible types are:

* transportOperatorST, performs surface tracking (ST) or ray tracing
* transportOperatorDT, performs Woodcock delta tracking (DT)
* transportOperatorHT, performs a hybrid between ST and DT

  - cutoff (*optional*, default = 0.9): cutoff between ST and DT. If, at the particle
    energy, the ratio between the local material cross section and the majorant cross
    section is larger than the cutoff, DT is used; otherwise ST is used.

Example: ::

      transportOperator { type transportOperatorHT; cutoff 0.85; }

Collision Operator
------------------

The **collision operator** process all collision types. It samples the colliding nuclide
and the reaction, and calculates all relevant by-products. In the input file, one must
include: ::

      collisionOperator { neutronCE { type <ceCollisionOperatorType>; *keywords* } }

if continuos energy nuclear data are used, or ::

      collisionOperator { neutronMG { type <ceCollisionOperatorType>; } }

if multi-group nuclear data are used. In a hybrid simulation, both ``neutronCE`` and
``neutronMG`` can be included.

The possible types to be used with **continuous energy** data are:

neutronCEstd
############

neutronCEstd, to perform analog collision processing

* minEnergy (*optional*, default = 1.0e-11): minimum energy cut-off. [MeV]
* maxEnergy (*optional*, default = 20.0): maximum energy cut-off. [MeV]
* energyThreshold (*optional*, default = 400): energy threshold for explicit treatment
  of target nuclide movement. Target movement is sampled if neutron energy E < kT ∗
  energyThreshold where kT is target material temperature in [MeV]. [-]
* massThreshold (*optional*, default = 1): mass threshold for explicit treatment of
  target nuclide movement. Target movement is sampled if target mass A < massThreshold. [Mn]

Example: ::

      collisionOperator { neutronCE { type neutronCEstd; minEnergy 1.0e-12; maxEnergy 30.0;
      energyThreshold 200; massThreshold 2; } }

neutronCEimp
############

neutronCEimp, to perform implicit collision processing

* minEnergy (*optional*, default = 1.0e-11): minimum energy cut-off. [MeV]
* maxEnergy (*optional*, default = 20.0): maximum energy cut-off. [MeV]
* energyThreshold (*optional*, default = 400): energy threshold for explicit treatment
  of target nuclide movement. Target movement is sampled if neutron energy E < kT ∗
  energyThreshold where kT is target material temperature in [MeV]. [-]
* massThreshold (*optional*, default = 1): mass threshold for explicit treatment
  of target nuclide movement. Target movement is sampled if target mass A <
  massThreshold. [Mn]
* splitting (*optional*, default = 0): 1 for true; 0 for false; enables splitting
  for particles above a certain weight
* roulette (*optional*, default = 0): 1 for true; 0 for false; enables rouletting
  of particles below a certain weight
* minWgt (*optional*, default = 0.25): minimum particle weight for rouletting
* maxWgt (*optional*, default = 1.25): maximum particle weight for splitting
* avgWgt (*optional*, default = 0.5): weight of a particle on surviving rouletting
* impAbs (*optional*, default = 0): 1 for true; 0 for false; enables implicit capture
* impGen (*optional*, default = 1): 1 for true; 0 for false; enables implicit fission
  sites generation
* weightWindows (*optional*, default = 0): 1 for true; 0 for false; enables the use of
  weight windows
* UFS (*optional*, default = 0): 1 for true; 0 for false; enables the use of uniform
  fission sites

Example: ::

      collisionOperator { neutronCE { type neutronCEimp; minEnergy 1.0e-12; maxEnergy 30.0;
      impAbs 1; roulette 1; splitting 1; impGen 1; maxWgt 2.0; minWgt 0.1; UFS 1; } }

The possible types to be used with **multi-group** data are:

neutronMGstd
############

neutronMGstd, to perform analog collision processing

Example: ::

      collisionOperator { neutronMG { type neutronMGstd; } }

Weight Windows
--------------

Weight windows can be used if, inside the collision operator ``CEneutronimp``, the
keyword ``weightWindows`` is set to 1. Then, in the input file, one needs to add: ::

        varianceReduction { type weightWindowsField; file <pathToWeightWindowsFile>; }

The file that contains **weight windows** has to include:

* map: map as defined for the tallies
* wLower: array with the lower weight windows weights, where the order of the values
  in the array must correspond to the order of the bins in the map
* wUpper: array with the upper weight windows weights, where the order of the values
  in the array must correspond to the order of the bins in the map
* constSurvival: multiplication constant. Multiplied by the lower weights, gives the
  survival weight for Russian roulette

Example: ::

      map  { type multiMap; maps (mapx mapy);
      mapx { type spaceMap;  axis x;  grid unstruct;  bins (0.0 1.0 2.0); }
      mapy { type spaceMap;  axis y;  grid unstruct;  bins (0.0 5.0 10.0 15.0); } }
      constSurvival 2.0;
      wLower (0.5 0.1 0.2 0.1 0.5 0.5);
      wUpper (2.0 1.2 1.5 1.1 2.0 4.0);

Uniform Fission Sites
---------------------

Weight windows can be used if, inside the collision operator ``CEneutronimp``, the
keyword ``UFS`` is set to 1. Then, in the input file, one needs to add: ::

      uniformFissionSites { type uniFissSitesField; map { <Map definition> } *keywords* }

In the input above, ``map`` is the geometrical map used for UFS. The map has to contain
fissile material for the method to make sense. Other keywords are:

* uniformVolMap (*optional*, default = 1): 1 for true; 0 for false; flag that states
  whether the bins of the map contain equal volumes of fissile material or not
* popVolumes (*optional*, default = 1.0e7): if ``uniformVolMap`` is false, a Monte Carlo
  calculation is run to estimate the fissile material volumes in each map bin. This entry
  correspond to the number of points sampled in the geometry for the volume calculation.
  Note that this volume calculation is done only once during initialisation

Example: ::

      uniformFissionSites { type uniFissSitesField; uniformVolMap 0; popVolumes 1.0e8;
      map { <Map definition> }
      }

Geometry
--------

A detailed description about the geometry modelling adopted in SCONE can be found at
:ref:`Geometry <Geometry>`. In an input file, one has to include: ::

      geometry  { type <geometryType>; boundary (a b c d e f); graph { type <graphType>; }
      surfaces  { <Surfaces definition> }
      cells     { <Cells definition> }
      universes { <Universes definition> }
      }

At the moment, the only **geometry** type available is ``geometryStd``. As for the boundary
six integers have to be inputted. These correspond to the boundary conditions at boundaries
(-x +x -y +y -z +z). The possibilities are:

* vacuum, or black: input 0
* reflective: input 1
* periodic: input 2

.. note::
    Strictly speaking it is up to a particular boundary surface to interpret how the values
    in the boundary condition sequence are interpreted. For all cube-like surfaces the rule
    above holds, but for more exotic boundaries (e.g., hexagons) it is worth double checking
    the documentation comment of the particular surface in the source code.

.. note::
   Curved surfaces only allow for vacuum boundaries.

The **graph** definition allows two options:

* shrunk: each local (material) cell has the same uniqueID in all universe instances
* extended: every local (material) cell has its own uniqueID in all universe instances

Hence, an example of a geometry input could look like: ::

      geometry  { type geometryStd; boundary (1 1 1 1 0 0); graph { type shrunk; }
      surfaces  { <Surfaces definition> }
      cells     { <Cells definition> }
      universes { <Universes definition> }
      }

For more details about the graph-like structure of the nested geometry see the relevant
:ref:`section <DAG_GEOM>`.

Surfaces
########

To define one or multiple **surfaces**, the necessary entries are: ::

      surfaces {
      <name1> { id <idNumber1>; type <surfaceType>; *keywords* }
      <name2> { id <idNumber2>; type <surfaceType>; *keywords* }
      ...
      <nameN> { id <idNumberN>; type <surfaceType>; *keywords* }
      }

Here, the ``name`` can be anything at the discretion of the user, as long as it doesn't
contain spaces. The ``idNumber`` can be any integer; attention must be paid that all
``idNumbers`` are unique.

Several ``surfaceTypes`` are possible:

* box: axis aligned box

  - origin: (x y z) vector with the origin position. [cm]
  - halfwidth: (x y z) vector with the halfwidth of each side. [cm]

Example: ::

      surf1 { id 92; type box; origin (0.0 0.0 9.0); halfwidth (1.0 2.0 0.3); }

* squareCylinder: infinitely long square cylinder aligned with x, y or z axis. The
input type has to be ``xSquareCylinder``, ``ySquareCylinder`` or ``zSquareCylinder``

  - origin: (x y z) vector with the origin position; the entry corresponding to
    the cylinder axis is ignored. [cm]
  - halfwidth: (x y z) vector with the halfwidth of each side; the entry
    corresponding to the cylinder axis is ignored. [cm]

Example: ::

      surf2 { id 25; type ySquareCylinder; origin (3.0 0.0 9.0); halfwidth (4.4 0.0 0.1); }

* truncCylinder: finite length cylinder aligned with x, y or z axis. The input
  type has to be ``xTruncCylinder``, ``yTruncCylinder`` or ``zTruncCylinder``

  - origin: (x y z) vector with the origin position. [cm]
  - halfwidth: axial halfwidth. [cm]
  - radius: cylinder radius. [cm]

Example: ::

      surf3 { id 3; type zTruncCylinder; origin (3.0 2.1 5.0); halfwidth 20.0;
      radius 1.6; }

* aPlane: plane with normal along x, y or z. The input type has to be ``xPlane``,
  ``yPlane`` or ``zPlane``

  - a0: position of the plane on the axis. The input type has to be ``x0``, ``y0``
    or ``z0``. [cm]

Example: ::

      surf4 { id 8; type xPlane; x0 4.0; }

* plane: generic plane (F(r) = c1 * x + c2 * y + c3 * z - c4)

  - coeffs: (c1 c2 c3 c4) vector with coefficients

Example: ::

      surf5 { id 55; type plane; coeffs (8.6 3.0 66.0 1.5); }

* cylinder: infinitely long cylinder aligned with x, y or z axis. The input type
  has to be ``xCylinder``, ``yCylinder`` or ``zCylinder``

  - origin: (x y z) vector with the origin position; the entry corresponding to
    the cylinder axis is ignored. [cm]
  - radius: cylinder radius. [cm]

Example: ::

      billy { id 92; type xCylinder; origin (0.0 0.0 9.0); radius 4.8; }

* sphere

  - origin: (x y z) vector with the origin position. [cm]
  - radius: sphere radius. [cm]

Example: ::

      surf6 { id 234; type sphere; origin (5.0 86.0 19.4); radius 18.3; }

Cells
#####

Similarly to the surfaces, the **cells** in the geometry can be defined as: ::

      cells {
      <name1> { id <idNumber1>; type <cellType>; surfaces (<surfaces>); filltype <fillType>; *keywords* }
      <name2> { id <idNumber2>; type <cellType>; surfaces (<surfaces>); filltype <fillType>; *keywords* }
      ...
      <nameN> { id <idNumberN>; type <cellType>; surfaces (<surfaces>); filltype <fillType>; *keywords* }
      }

At the moment, in SCONE, the only ``cellType`` available is ``simpleCell``.
In the surface definition, one should include the indexes of the corresponding
surfaces with no sign to indicate a positive half-space, or minus sign to indicate
a negative half-space. The space in between cells corresponds to an intersection.

The possible ``fillTypes`` are:

* mat: if the cells is filled with a homogeneous material

  - material: takes as an input the material name

Example: ::

      cell1 { id 1; type simpleCell; surfaces (1 -6 90); filltype mat; material fuel; }

* uni: if the cell is filled with a universe

  - universe: takes as an input the universe ``id``

Example: ::

      cellX { id 5; type simpleCell; surfaces (2 -3); filltype uni; universe 6; }

* outside: if the cell is outside of the geometry

Example: ::

      cellixx { id 55; type simpleCell; surfaces (-10); filltype outside; }

Universes
#########

Similarly to the surfaces and cells, the **universes** in the geometry can be defined as: ::

      universes {
      <name1> { id <idNumber1>; type <universeType>; *keywords* }
      <name2> { id <idNumber2>; type <universeType>; *keywords* }
      ...
      <nameN> { id <idNumberN>; type <universeType>; *keywords* }
      }

Several ``universeTypes`` are possible:

* cellUniverse, composed of the union of different cells. Note that overlaps are
  forbidden, but there is no check to find overlaps

  - cells: array containing the ``cellIds`` as used in the cell definition
  - origin (*optional*, default = (0.0 0.0 0.0)): (x y z) array with the origin
    of the universe. [cm]
  - rotation (*optional*, default = (0.0 0.0 0.0)): (x y z) array with the
    rotation angles in degrees applied to the universe. [°]

.. note::
   When creating a ``cellUniverse`` a user needs to take care to avoid leaving
   any 'unspecified' regions (sets in space which do not belong to any cell).
   If these are reachable by a particle (e.g., are not covered by any higher
   level universe) they will cause a calculation to crash.

Example: ::

      uni3 { id 3; type cellUniverse; cells (1 2 55); origin (1.0 0.0 0.0); rotation (0.0 90.0 180.0); }

* pinUniverse, composed of infinite co-centred cylinders

  - radii: array containing the radii of the co-centred cylinders. There
    must be an entry equal to 0.0, which corresponds to the outermost
    layer, which is infinite. [cm]
  - fills: array containing the names or ids of what is inside each cylindrical
    shell. The order of the fills must correspond to the order of the corresponding
    radii. An entry can be a material name, the keyword ``void``, or a   ``u<id>``,
    where ``id`` is the id of a defined universe
  - origin (*optional*, default = (0.0 0.0 0.0)): (x y z) array with the
    origin of the universe. [cm]
  - rotation (*optional*, default = (0.0 0.0 0.0)): (x y z) array with the
    rotation angles in degrees applied to the universe. [°]

Example: ::

      uni3 { id 3; type pinlUniverse; radii (0.2 1.0 1.1 1.3 0.0); fills (u<1> fuel void clad coolant); }

* latUniverse, cartesian lattice of constant pitch

  - shape: (x y z) array of integers, stating the numbers of x, y and z
    elements of the lattice. For a 2D lattice, one of the entries has to be 0
  - pitch: (x y z) array with the x, y and z lattice pitches. In a 2D lattice,
    the value entered in the third dimension is not used. [cm]
  - padmat: material name or universe index (u<id>) that fills the possible
    extra space between the lattice and its bounding surface. Also the keyword
    ``void`` is allowed
  - map: map that includes the universe ids of the elements of the lattice.
    The order is: increasing x, increasing y and then increasing z
  - origin (*optional*, default = (0.0 0.0 0.0)): (x y z) array with the
    origin of the universe. [cm]
  - rotation (*optional*, default = (0.0 0.0 0.0)): (x y z) array with the
    rotation angles in degrees applied to the universe. [°]

Example: ::

      uni_lattice { id 10; type latUniverse; shape (3 2 2); pitch (1.0 1.0 1.5); padMat u<3>; map (
      1 2 3 // x: 1-3, y: 2, z: 2
      4 5 6 // x: 1-3, y: 1, z: 2
      7 8 9 // x: 1-3, y: 2, z: 1
      10 11 12 ) } // x: 1-3, y: 1, z: 1

.. note::
   The order of the elements in the lattice is different from other MC codes, e.g.,
   Serpent. The lattice is written in the style *WYSIWYG*: What You See Is What You Get.

* rootUniverse: top level universe of geometry

  - border: id of the boundary surface for the whole geometry
  - fill: inside filling, as a material name or a universe (u<id>)

Example: ::

      root { id 1000; type rootUniverse; border 10; fill u<1>; }

Visualiser
----------

To **plot** a geometry, the keyword ``viz`` must be present in the input file: ::

      viz {
      <name1> { type <vizType>; *keywords* }
      <name2> { type <vizType>; *keywords* }
      }

The possible types of files that the geometry is plotted in are:

vtk
###

* corner: (x y z) array with the corner of the geometry [cm]
* width: (x y z) array with the width of the mesh in each direction [cm]
* vox: (x y z) array with the number of voxels requested in each direction
* what (*optional*, default = material): defines what is highlighted in the
  plot; options are ``material`` and ``cellID``

Example: ::

      plotVTK { type vtk; corner (10.0 6.0 2.0); width (20.0 12.0 4.0); vox (4000 120 400); what cellID; }

bmp
###

* centre: (x y z) array with the coordinates of the center of the plot [cm]
* axis: ``x``, ``y`` or ``z``, it's the axis normal to the 2D plot
* width (*optional*, default = whole geometry): (y z), (x z) or (x y) array
  with the width of the geometry plotted in each direction [cm]
* res: (y z), (x z) or (x y) array with the resolution of the mesh in each direction
* output: name of the output file, with extension ``.bmp``
* what (*optional*, default = material): defines what is highlighted in the
  plot; options are ``material`` and ``cellID``

Example: ::

      plotBMP { type bmp; axis z; width (50 10); res (1000 200); output geomZ; what material; }

.. note::
   SCONE can be run to visualise geometry without actually doing transport, by
   including ``--plot`` when running the application. In this case the visualiser
   has to be included in the file.

Nuclear Data
------------

SCONE can be used with both continuous energy data and multi-group data. The type
of data used must be specified in the ``physicsPackage`` options, as well as in the
``collisionOperator`` options. As for **nuclear data**, the input files has to look like: ::

      nuclearData {
      handles { <Nuclear data handles definition> }
      materials { <Materials definition> }
      }

The **handles** definition is structured as the following: ::

      handles {
      <handleName1> { type <databaseType>; *keywords* }
      <handleName2> { type <databaseType>; *keywords* }
      }

The name of a handle has to be the same as defined in a ``physicsPackage`` under the
keyword ``XSdata``.

Otherwise, the possible **nuclear database** types allowed are:

aceNeutronDatabase
##################

aceNeutronDatabase, used for continuous energy data. In this case, the data is read
from ACE files.

* aceLibrary: includes the path to the *.aceXS* file, which includes the paths to
  the ACE files
* ures (*optional*, default = 0): 1 for true; 0 for false; activates the unresolved
  resonance probability tables treatment

Example: ::

      ceData { type aceNuclearDatabase; aceLibrary ./myFolder/ACElib/JEF311.aceXS; ures 1; }

baseMgNeutronDatabase
#####################

baseMgNeutronDatabase, used for multi-group data. In this case, the data is read
from files provided by the user.

* PN: includes a flag for anisotropy treatment. Could be ``P0`` or ``P1``

Example: ::

      mgData { type baseMgNeutronDatabase; PN P1; }

The *materials* definition is structured as: ::

      materials {
      <materialName1> { temp <temp1>;
      composition { <Composition definition> }
      *keywords* }
      <materialName2> { temp <temp2>;
      composition { <Composition definition> }
      *keywords* }
      }

In this case, ``materialName`` can be any name chosen by the user; ``temp`` is the
material temperature in [K].

The ``composition`` dictionary must always be included, but it can be empty in
multi-group simulations. In continuous energy simulations, it should include a
list of the ZAIDs of all the nuclides that compose that material, and the respective
atomic densities in [atoms/cm/barn]. The ZAIDs are normally in the form ``ZZAAA.TT``,
or ``ZAAA.TT`` for nuclides with Z<10. The code ``TT`` indicates the temperature used
in the nuclear data evaluation, and the options are 03, 06, 09, 12 and 15,
corresponding to temperatures of 300K, 600K, 900K, 1200K and 1500K.

Other options are:

* moder: dictionary that includes information on thermal scattering data. It has to
  include a list of ZAIDs for which S(a,b) has to be used, and the name of the file
  that contains the data. The file has to be included in the list of files in the *.aceXS*
  input file. Note that this input is ignored if the nuclide or nuclides listed are not
  included in the material. Only needed for continuous energy simulations.

* xsFile: needed for multi-group simulations. Must contain the path to the file where
  the multi-group cross sections are stored.

Example 1: ::

      materials {
      fuel { temp 273;
      composition {
      92238.03   0.021;
      92235.03   0.004;
      8016.03    0.018535464; }
      }
      water { temp 273;
      composition {
      1001.03   0.0222222;
      8016.03   0.00535; }
      moder { 1001.03 h-h2o.42; }
      }
      }

Example 2: ::

      materials {
      fuel { temp 573;
      composition { }
      xsFile ./xss/fuel.txt
      }
      }

Multi-group cross sections
--------------------------

In the case of a multi-group calculation, **multi-group cross sections** must be
provided by the user. These are in separate files compared to the input file. The
structure of such cross section files is the following: they must include

* numberOfGroups: number of energy groups used (=N)
* capture: vector of size N with the material-wise macroscopic capture cross section.
  The order of the elements corresponds to groups from fast (group 1) to thermal
  (group N)
* fission (*optional*): vector of size N with the material-wise macroscopic fission
  cross section. The order of the elements corresponds to groups from fast (group 1)
  to thermal (group N). Must be included only if the materials is fissile
* nu (*optional*): vector of size N with the material-wise macroscopic neutron
  production nu-bar. The order of the elements corresponds to groups from
  fast (group 1) to thermal (group N). Must be included only if the materials
  is fissile
* chi (*optional*): vector of size N with the material-wise fission spectrum. The order
  of the elements corresponds to groups from fast (group 1) to thermal (group N).
  Must be included only if the materials is fissile
* P0: P0 scattering matrix, of size NxN. In the case of a 3x3 matrix, the elements are
  ordered as: ::

      1 -> 1   1 -> 2   1 -> 3
      2 -> 1   2 -> 2   2 -> 3
      3 -> 1   3 -> 2   3 -> 3

* scatteringMultiplicity: P0 scattering multiplicity matrix, of size NxN. Contains
  multiplicative elements that will be multiplied to the P0 matrix elements for scattering
  production cross section, hence all elements must be >= 1.0
* P1 (*optional*): necessary only if ``P1`` is defined in the ``baseMgNeutronDatabase``
  entry ``PN``. It contains the P1 scattering matrix, of size NxN

An example file is: ::

      numberOfGroups 2;
      capture (0.0010046 0.025788);
      fission (0.0010484 0.050632);
      nu      (2.5 2.5);
      chi     (1.0 0.0);
      scatteringMultiplicity (
      1.0 1.0
      1.0 1.0  );
      P0 (
      0.62568 0.029227
      0.0     2.443830
      );
      P1 (
      0.27459 0.0075737
      0.0     0.83318
      );

Tallies
-------

As mentioned previously, one might have to include the keywords ``inactiveTally`` and
``activeTally`` in the input file (in the case of ``eigenPhysicsPackage``), or just
``tally`` (in the case of ``fixedSourcePhysicsPackage``). Either way, the **tally**
definition is the same for all cases: ::

      tally {
      *keywords*
      <resName1> { type <clerkType1>; response (<responseName>); <responseName> { type <responseType>; *keywords* } *keywords* }
      <resName2> { type <clerkType2>; *keywords* }
      ...
      <resNameN> { type <clerkTypeN>; }
      }

In this case, ``resName`` can be any name chosen by the user, and it is what will be
reported in the output file.

Tally Clerks
############

The **tally clerks** determine which kind of estimator will be used. The options are:

* collisionClerk, for a collision estimator of flux and reaction rates

  - response: defines which response function has to be used for this tally. Note
    that more than one response can be defined per each tally
  - map (*optional*): contains a dictionary with the ``tallyMap`` definition,
    that defines the domains of integration of each tally
  - filter (*optional*): can filter out particles with certain properties,
    preventing them from scoring results

* trackClerk

  - response: defines which response function has to be used for this tally.
    Note that more than one response can be defined per each tally
  - map (*optional*): contains a dictionary with the ``tallyMap`` definition,
    that defines the domains of integration of each tally
  - filter (*optional*): can filter out particles with certain properties,
    preventing them from scoring results

Example: ::

      tally {
      collision_estimator { type collisionClerk; response (<responseName>); <responseName> { type <responseType>; *keywords* }
      map { <Map definition> }
      filter { <Filter definition> }
      }
      track_estimator { type trackClerk; response (<responseName1> <responseName2>);
      <responseName1> { type <responseType>; *keywords* }
      <responseName2> { type <responseType>; *keywords* }
      }
      }

* keffAnalogClerk, analog k_eff estimator
* keffImplicitClerk, implicit k_eff estimator

Example: ::

      tally {
      k_eff1 { type keffAnalogClerk; }
      k_eff2 { type keffImplicitClerk; }
      }

* centreOfMassClerk, geometrical 3D center of mass estimator

  - cycles: number of cycles for which to track center of mass

Example: ::

      tally {
      com { type comClerk; cycles 200; }
      }

* collisionProbabilityClerk, tallies a collision probability matrix

  - map: contains a dictionary with the ``tallyMap`` definition, that defines
    the bins of the matrix

Example: ::

      tally {
      collisionProb { type collisionProbabilityClerk; map { <Map definition> } }
      }

* dancoffBellClerk, calculates a single-term rational approximation for a lattice

  - fuelMat: list of fuel material names
  - modMat: list of moderator material names
  - Elow (*optional*, default = 0.0): bottom energy boundary; [MeV]
  - Etop (*optional*, default = 20.0): top energy boundary; [MeV]

Example: ::

      tally {
      dancoff_bell_factors { type dancoffBellClerk; fuelMat (fuel1 fuel2 fuel_Gd); modMat (water); Elow 0.06; Etop 10.0; }
      }

* mgXsClerk, calculates multi-group cross sections via a collision estimator
  of reaction rates and analog tallies of fission spectrum and scattering events
  ingoing and outgoing energies and multiplicity

  - energyMap (*optional*, default = 1 group): definition of the energy group
    structure to be used
  - spaceMap (*optional*, default = whole geometry): definition of a spatial
    tally map
  - PN (*optional*, default = 0): 1 for true; 0 for false; flag that indicates
    whether to calculate scattering matrices only up to P1 (``PN 0``) or P7 (``PN 1``)

Example: ::

      tally {
      MGxss { type mgXsClerk;
      energyMap { <Map definition> }
      spaceMap { <Map definition> }
      PN 1; }
      }

* shannonEntropyClerk, implicit Shannon entropy estimator

  - map: contains a dictionary with the ``tallyMap`` definition, that defines
    the (spatial) discretisation used to score the entropy
  - cycles: number of cycles to tally the entropy for

Example: ::

      tally {
      shannon_entropy { type shannonEntropyClerk;
      map { <Map definition> }
      cycles 200; }
      }

* simpleFMClerk, 1D fission matrix collision estimator

  - map: contains a dictionary with the ``tallyMap`` definition, that defines
    the bins of the matrix

Example: ::

      tally {
      fissionMat { type simpleFMClerk; map { <Map definition> } }
      }

Tally Responses
###############

Certain tally clerks, like the ``collisionClerk`` and ``trackClerk``, require
a **response function**. The different types of responses could be:

* fluxResponse: used to calculate the flux, i.e., the response function is 1.0

Example: ::

      tally {
      collision_estimator { type collisionClerk; response (flux); flux { type fluxResponse; } }
      }

* macroResponse: used to score macroscopic reaction rates

  - MT: MT number of the desired reaction. The options are: -1 total, -2 capture,
    -6 fission, -7 nu*fission, -21 absorption

Example: ::

      tally {
      collision_estimator { type collisionClerk; response (total fission);
      total { type macroResponse; MT -1; }
      fission { type macroResponse; MT -6; } }
      }

* microResponse: used to score microscopic reaction rates

  - MT: MT number of the desired reaction. The options are: 1 total, 2 elastic
    scattering, 18 fission, 27 absorption, 102 capture
  - material: material name where to score the reaction. The material must be
    defined to include only one nuclide; its density could be anything, it doesn't
    affect the result

Example: ::

      tally {
      collision_estimator { type collisionClerk; response (elScatter capture);
      elScatter { type microResponse; MT 2; material water; }
      capture { type microResponse; MT 102; material fuel; }
      }
      }

* weightResponse: response for scoring particle weights

  - moment (*optional*, default = 1): moment of the weight scored

Example: ::

      tally {
      collision_estimator { type collisionClerk; response (weight0 weight1 weight2);
      weight0 { type weightResponse; moment 0; }
      weight1 { type weightResponse; moment 1; }
      weight2 { type weightResponse; moment 2; }
      }
      }

.. note::
   To calculate the average weight, one should divide weight moment 1 (weight1)
    by weight moment 0 (weight0). To calculate the variance of the weights, the
    tally results have to be post-processed as: var = weight2/weight0 - (weight1/weight0)^2

Tally Maps
##########

The different types of **tally maps** are:

* cellMap (1D map), cell-wise map

  - cells: list of ids of the cells to be used an map bins
  - undefBin (*optional*, default = false): 'yes','y','true','TRUE','T' for true;
    'no', 'n', 'false', 'FALSE', 'F' for false; flag that indicates whether all
    the cells not listed in ``cells`` should constitute a map bin or not

Example: ::

      map { type cellMap; cells (1 5 3 2 4 100); undefBin T; }

* energyMap (1D map), defines an energy group structure

  - grid: ``log`` for logarithmically spaced bins or ``lin`` for linearly spaced bins

    + min: bottom energy [MeV]
    + max: top energy [MeV]
    + N: number of bins

  - grid: ``unstruct`` for unstructured grids, to be manually defined

    + bins: array with the explicit definition of the energy bin boundaries to be used

  - grid: ``predef``

    + name: name of the predefined group structure. Options are: ``wims69``,
      ``wims172``, ``casmo40``, ``casmo23``, ``casmo12``, ``casmo7``, ``vitaminj``

Examples: ::

      map1 { type energyMap; grid log; min 1.0e-11; max 20.0; N 300; }
      map2 { type energyMap; grid lin; min 1.0; max 20.0; N 100; }
      map3 { type energyMap; bins (1.0E-9 1.0E-8 0.6E-6 0.3 20.0); }
      map4 { type energyMap; name casmo12; }

* homogMatMap (1D map), divides based on the material a particle is in with the
  possibility of grouping some materials together

  - bins: list of names of the material bins, that can contain one or more
    materials; this is followed by all the bin names as key, and the material
    names included in the bin as an entry
  - undefBin (*optional*, default = false): 'yes','y','true','TRUE','T' for true;
    'no', 'n', 'false', 'FALSE', 'F' for false; flag that indicates whether all
    the materials not included in any bin should constitute a map bin or not

Example: ::

      map { type homogMatMap; bins (bin1 bin2 bin3);
      bin1 (mat1 mat2 mat3);
      bin2 (fuel1 fuel3 uo2);
      bin3 (water);
      undefBin T;
      }

* materialMap (1D map), material-wise map

  - materials: list of material names to be used as map bins
  - undefBin (*optional*, default = false): 'yes','y','true','TRUE','T' for true;
    'no', 'n', 'false', 'FALSE', 'F' for false; flag that indicates whether all
    the materials not included should constitute a map bin or not

Example: ::

      map { type materialMap; materials (fuel water cladding reflector fuelGd); undefBin T; }

* multiMap, ensemble of multiple 1D maps

  - maps: list of the names of the maps that will compose the ``multiMap``. This
    is followed by dictionaries that define the requested maps

Example: ::

      map { type multiMap; maps (map1 map2 map10);
      map1 { <1D map definition> }
      map2 { <1D map definition> }
      map10 { <1D map definition> }
      }

* spaceMap (1D map), geometric cartesian map

  - axis: ``x``, ``y`` or ``z``

  - grid: ``lin`` for linearly spaced bins

    + min: bottom coordinate [cm]
    + max: top coordinate [cm]
    + N: number of bins

  - grid: ``unstruct`` for unstructured grids, to be manually defined

    + bins: array with the explicit definition of the bin boundaries to be used

Examples: ::

      map1 { type spaceMap; axis x; grid lin; min -50.0; max 50.0; N 100; }
      map2 { type spaceMap; axis z; grid unstruct; bins (0.0 0.2 0.3 0.5 0.7 0.8 1.0); }

* sphericalMap, geometric spherical map

  - origin (*optional*, default = (0.0 0.0 .0.)): (x y z) vector with the origin
    of the spherical map

  - grid: ``lin`` for linearly spaced bins or ``equivolume`` for spherical shells

    + Rmin (*optional*, default = 0.0): minimum radius [cm]
    + Rmax: maximum radius [cm]
    + N: number of radial bins

  - grid: ``unstruct`` for unstructured grids, to be manually defined

    + bins: array with the explicit definition of the spherical bin boundaries
      to be used

Examples: ::

      map1 { type sphericalMap; origin (2.0 1.0 0.0); grid lin; Rmin 3.0; Rmax 10.0; N 14; }
      map2 { type sphericalMap; grid equivolume; Rmax 20.0; N 10; }
      map3 { type sphericalMap; grid unstruct; bins (1.0 2.0 2.5 3.0 5.0); }

* cylindricalMap, geometric cylindrical map; other than the radial discretisation,
  one could add axial and azimuthal discretisation

  - orientation (*optional*, default = ``z``): ``x``, ``y`` or ``z``, axial direction
  - origin (*optional*, default = (0.0 0.0)): (y z), (x z) or (x y) vector with
    the origin of the cylindrical map
  - rGrid: ``lin`` for linearly spaced bins or ``equivolume`` for cylindrical shells

    + Rmin (*optional*, default = 0.0): minimum radius [cm]
    + Rmax: maximum radius [cm]
    + rN: number of radial bins

  - rGrid: ``unstruct`` for unstructured grids, to be manually defined

    + bins: array with the explicit definition of the cylindrical radial bin
      boundaries to be used

  - axGrid (*optional*, default = 1 bin): ``lin`` for linearly spaced axial bins

    + axMin: minimum axial coordinate [cm]
    + axMax: maximum axial coordinate [cm]
    + axN: number of axial bins

  - azimuthalN (*optional*, default = 1 bin): number of angular azimuthal bins

Example: ::

      map1 { type cylindricalMap; orientation y; origin (7.0 0.0); rGrid lin; Rmax 5.0; rN 10; }
      map2 { type cylindricalMap; rGrid unstruct; bins (2.0 3.0 4.5 5.0); axGrid lin; axMin 0.0; axMax 6.0 axN 24; azimuthalN 8; }

* weightMap (1D map), divides weight into number of discrete bins

  - grid: ``log`` for logarithmically spaced bins or ``lin`` for linearly spaced bins

    + min: bottom weight
    + max: top weight
    + N: number of bins

  - grid: ``unstruct`` for unstructured grids, to be manually defined

    + bins: array with the explicit definition of the weight bin boundaries to be used

Examples: ::

      map1 { type weightMap; grid log; min 1.0e-3; max 100.0; N 100; }
      map2 { type weightMap; grid lin; min 0.1; max 2.0; N 20; }
      map3 { type weightMap; bins (0.0 0.2 0.4 0.6 0.8 1.0 2.0 5.0 10.0); }

Tally Filters
#############

Another option that can be included in the tallies is **tally filters**. These
allow to filter out certain types of particles when scoring results. For now,
the only type of filter existing is:

* energyFilter, to stop particles within a certain energy range from contributing
  to a certain tally

  - Emin (for continuous energy particles): minimum energy [MeV]
  - Emax (for continuous energy particles): maximum energy [MeV]
  - Gtop (for multi-group particles): top energy group
  - Glow (for multi-group particles): bottom energy group

Example: ::

      CEfilter { type energyFilter; Emin 10.0; Emax 20.0; }
      MGfilter { type energyFilter; Gtop 1; Glow 5; }

Other options
#############

Other keywords, such as for results **normalisation**, that could be included are:

* norm: its entry is the name of the tally, ``resName``, to be used as a normalisation
  criterion. If the tally has multiple bins, (e.g. has a map), the bin with index 1
  will be used for normalisation
* normVal: value to normalise the tally ``resName`` to
* display: its entry is the name of the tally, ``resName``, which will be displayed
  each cycle. Only the tally clerks ``keffAnalogClerk`` and ``keffImplicitClerk``
  support display at the moment
* batchSize (*optional*, default = 1): the number of cycles that constitute a single
  batch for the purpose of statistical estimation. For example, a value of 5 means
  that a single estimate is obtained from a score accumulated over 5 cycles

Example: ::

      tally  {
      display (k-eff);
      norm fissRate;
      normVal 100.0;
      k-eff { type keffAnalogClerk;}
      fissRate { type collisionClerk; response (fission); fission {type macroResponse; MT -6;} }
      }
