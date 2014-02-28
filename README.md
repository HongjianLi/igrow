igrow
=====

igrow is a multithreaded structure-based [drug design] tool for computational drug discovery. It is inspired by [AutoGrow], and is hosted by GitHub at https://github.com/HongjianLi/igrow under [Apache License 2.0].


Features
--------

* igrow uses [idock] as backend docking engine.
* igrow implements branch exchange as its crossover operator.
* igrow digests ligands in pdbqt format, saving the effort of calling the prepare_ligand4 python script.
* igrow uses io service pool to reuse threads and maintain a high CPU utilization throughout the entire synthsizing procedure.
* igrow utilizes dynamic pointer vector to cache and sort ligands.
* igrow traces the sources of generated ligands and dumps the statistics in csv format so that users can easily get to know how the ligands are synthesized from the initial elite ligands.


Supported operating systems and compilers
-----------------------------------------

* Arch Linux 3.12.9 x86_64 and CLANG 3.4
* Windows 7 SP1 x64 and Visual Studio 2013 Update 1


Compilation
-----------

igrow depends on [Boost C++ Libraries]. Boost 1.55.0 is supported. The must-be-built libraries required by igrow are `System`, `Filesystem` and `Program Options`. An unofficial header-only library, Boost.Process, is also required by igrow. The file `process.zip` must be extracted to the Boost distribution tree in order to pass compilation.

### Compilation on Linux

The Makefile uses clang as the default compiler. To compile, simply run

    make

One may modify the Makefile to use a different compiler or different compilation options.

The generated objects will be placed in the `obj` folder, and the generated executable will be placed in the `bin` folder.

### Compilation on Windows

Visual Studio 2013 solution and project files are provided. To compile, simply run

    msbuild /t:Build /p:Configuration=Release

Or one may open `igrow.sln` in Visual Studio 2013 and do a full rebuild.

The generated objects will be placed in the `obj` folder, and the generated executable will be placed in the `bin` folder.


Usage
-----

First add igrow to your PATH environment variable.

To display a full list of available options, simply run the program without arguments

    igrow

The `examples` folder contains several use cases. For example, to grow TMC278 and dock the generated ligands against HIV-RT of PDB ID 2ZD1,

    cd examples/2ZD1

One can supply the options from command line arguments

    igrow --initial_generation_csv ../../../idock/examples/2ZD1/ZINC/log.csv --fragment_folder ../../fragments --idock_config idock.cfg

Or one can instruct igrow to load the options from a configuration file

    igrow --config igrow.cfg


Documentation Creation
----------------------

Documentations in both HTML and LaTeX formats can be esaily created by running [doxygen]

    doxygen igrow.dox

The created documents will be placed in `doc` folder. To compile LaTeX files into PDF, one must have `pdflatex` installed.

    cd doc/latex
    make

The generated PDF will be `refman.pdf`.


Change Log
----------

### 1.0 (under construction)

* Used idock as backend docking engine.
* Supported direct PDBQT manipulation without file format conversion.
* Used dynamic pointer vector to cache ligands.
* Supported dumping statistics and traceability of created ligands.
* Used docked atom coordinates to construct child ligands of the next generation.
* Parallelized crossover operations.
* Provided precompiled executables for 64-bit Linux and Windows.


Author
--------------

[Jacky Lee]


Logo
----

![igrow logo](https://github.com/HongjianLi/igrow/raw/master/logo.png)


[drug design]: http://en.wikipedia.org/wiki/Drug_design
[AutoGrow]: http://autogrow.ucsd.edu
[idock]: https://github.com/HongjianLi/idock
[Apache License 2.0]: http://www.apache.org/licenses/LICENSE-2.0.html
[C++11]: http://en.wikipedia.org/wiki/C++11
[Boost C++ Libraries]: http://www.boost.org
[doxygen]: http://www.doxygen.org
[Jacky Lee]: http://www.cse.cuhk.edu.hk/~hjli
