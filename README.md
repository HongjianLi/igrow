- About
  SmartGrow is a drug design project which runs under linux.

- Setup script files

  Download AutoDock Vina, http://vina.scripps.edu/
  Download MGLTools, http://mgltools.scripps.edu/downloads
  Install the mentioned programs

autovina.dat

  "Absolute path to Vina" --receptor "<RECEPTORFILENAME>" --ligand "<LIGANDFILENAME>" --center_x <CENTER_X> --center_y <CENTER_Y> --center_z <CENTER_Z> --size_x <SIZE_X> --size_y <SIZE_Y> --size_z <SIZE_Z> --out "<COMBINEDFILENAME>" --log "<LOGFILENAME>"

  Change "Absolute path to Vina" to the system path to Autodock Vina

prepare_ligand.dat

  "Absolute Path/MGLTools 1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py" -l "<PDBFILENAME>" -o "<PDBQTFILENAME>"

  Change the "Absolute path" section such that it points to "prepare_ligand4.py" of MGLTools

prepare_receptor.dat

  "Absolute Path/MGLTools 1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py" -r "<FILENAME>" -o "<FILENAME>.pdbqt"

  Change the "Absolute path" section such that it points to "prepare_receptor4.py" of MGLTools

- Setting up SmartGrow

  Extract the packed resources files for scripts and fragment libraries

  Under Linux

  Unpack the source code, go to the root director of the source code and type "make all" in the command console


- Executing SmartGrow
  To execute SmartGrow, assuming the executable is called SmartGrow
  Type the following the in command console

    ./SmartGrow -run_mode Execute -parm_file default.prm

  It specified the running mode to be "Execute" and utilised "default.prm" as the parameters file.

  To rerun the program, remove all folders and files under the "working root directory".

- Parameter file
  Example parameter file

  //DIRECTORIES

  // The generated drug candidates are placed under this path
  working root directory: C:/Documents and Settings/Administrator/My Documents/Visual Studio 2010/Projects/SmartGrow/Debug

  // The fragment library to added to the initial ligand
  fragments directory: C:/Documents and Settings/Administrator/My Documents/Visual Studio   2010/Projects/SmartGrow/fragment

  // Necessary scripts to call Autodock Vina, conversion and other tools
  scripts directory: C:/Documents and Settings/Administrator/My Documents/Visual Studio   2010/Projects/SmartGrow/scripts

  //INPUT FILES

  // The path of the initial ligand
  initial ligand: C:/Documents and Settings/Administrator/My Documents/Visual Studio     2010/Projects/SmartGrow/HIV protease/3KFNlig.pdbqt

  // The path of target protein to dock
  receptor: C:/Documents and Settings/Administrator/My Documents/Visual Studio 2010/Projects/SmartGrow/HIV protease/3KFNrec.pdbqt

  //AUTODOCK PARAMETERS

  // initial orientation of docking parameters
  autodock grid center: 18.490 3.665 1.864
  autodock box size: 12 14 12

  //EVOLUTION PARAMETERS

  // The number of elitist to carryover to the next generation
  number of carryovers: 12

  // The number of individuals produced by crossover
  number of children:  8

  // The number of individuals produced by mutation
  number of mutants: 40

  // The maximum number of atoms an individual can substain
  max number atoms: 70

  // The set of indices of hydrogen that will not be used to add a fragment, please keep default = -1
  indices of hydrogens that are not linkers: -1

  // Total number of generation to execute, usually a multiple of "dock frequency"
  number of generations: 24

  // The period on how frequent the docking program is executed, 3 means docking program is called every 3 generations
  dock frequency: 3

  //OPTIONAL PARAMETERS

  // Control whether the evaluation to consider the molecular weight of the individuals, default = false
  score ligand size: false

  // Add hydrogens to the initial ligand when applicable, default = true
  add hydrogen: true

  // Advanced joining method can be used during mutation, default = true
  advanced synthesis mode: true

  // Incorporate drug-like properties in the constraint, default = true
  lipinski rule: true

  // The threshold parameter to decide how similar two individuals are, default = 0.01
  molecule similarity: 0.01