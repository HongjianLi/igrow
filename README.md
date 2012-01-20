igrow
=====

igrow is a multithreaded structure-based [drug design] tool for computational drug discovery. It is inspired by [AutoGrow], and is hosted by GitHub at https://github.com/HongjianLi/igrow under [Apache License 2.0].


Features
--------

* igrow uses either [idock] or [AutoDock Vina] as backend docking engine.
* igrow supports more types of chemical synthesis such as halogen replacement and branch replacement in addition to hydrogen replacement.
* igrow digests ligands and fragments in pdbqt format, saving the effort of frequently calling the prepare_ligand4 python script.
* igrow invents its own thread pool in order to reuse threads and maintain a high CPU utilization throughout the entire synhsizing procedure. The thread pool parallelizes the creation of mutants and children in each generation.
* igrow utilizes flyweight pattern for caching fragments and dynamic pointer vector for caching and sorting ligands.
* igrow traces the sources of generated ligands and dumps the statistics in csv format so that users can easily get to know how the ligands are synthesized from the initial ligand and fragments.


Tested operating systems and compilers
-----------------------------------------

* Ubuntu 11.10 x86_64 and GCC 4.6.1
* Ubuntu 11.10 x86_64 and CLANG 2.9
* Ubuntu 11.10 x86_64 and Intel C++ Compiler 12.0.5.220
* Fedora 16 x86_64 and GCC 4.6.2
* Fedora 16 x86_64 and Intel C++ Compiler 12.1.2
* Arch Linux 3.2.1 x86_64 and GCC 4.6.2
* Arch Linux 3.2.1 x86_64 and CLANG 3.0
* Arch Linux 3.2.1 x86_64 and Intel C++ Compiler 12.1.2
* Windows 7 SP1 x64 and Windows SDK 7.1
* Windows 7 SP1 x64 and Visual Studio 2010 SP1
* Windows 7 SP1 x64 and Intel C++ Compiler 12.1.2


Compilation
-----------

igrow depends on [Boost C++ Libraries]. Boost 1.48.0 is tested. The must-be-built libraries required by igrow are `System`, `Thread`, `Filesystem`, and `Program Options`. Two unofficial and header-only libraries, Boost.Process and Boost.Atomic, are also required by igrow. The file `boost.process.tar.bz2` and `boost.atomic.tar.bz2` must be extracted to the Boost distribution tree in order to pass compilation.

### Compilation on Linux

The Makefile uses GCC as the default compiler. To compile, simply run

    make

CLANG is also supported.

    make TOOLSET=clang

Intel C++ Compiler is also supported.

    make TOOLSET=intel

One may modify the Makefile to use a different compiler or different compilation options.

The generated objects will be placed in the `obj` folder, and the generated executable will be placed in the `bin` folder.

### Compilation on Windows

Visual Studio 2010 solution and project files are provided in the `msvc` folder. The project file uses Windows 7.1 SDK as platform toolset by default. One may revert it to vc100. To compile, simply run

    msbuild /t:Build /p:Configuration=Release

Or one may open `igrow.sln` in Visual Studio 2010 and do a full rebuild.

The generated objects will be placed in the `obj` folder, and the generated executable will be placed in the `bin` folder.


Usage
-----

First add igrow to your PATH environment variable.

To display a full list of available options, simply run the program without arguments

    igrow

The `examples` folder contains several use cases. For example, to grow TMC278 and dock the generated ligands against HIV-RT of PDB ID 2ZD1,

    cd examples/2ZD1

One can supply the options from command line arguments

    igrow --fragment_folder ../../fragments --initial_ligand ../../fragments/CH4.pdbqt --docking_program idock --docking_config docking.cfg

Or one can instruct igrow to load the options from a configuration file

    igrow --config igrow.cfg


Documentation Creation
----------------------

Documentations in both HTML and LaTeX formats can be esaily created by running [doxygen]

    doxygen doxygen

The created documents will be placed in `doc` folder. To compile LaTeX files into PDF, one must have `pdflatex` installed.

    cd doc/latex
    make

The generated PDF will be `refman.pdf`.


Change Log
----------

### 1.0

* Supported using either idock or AutoDock Vina as backend docking engine.
* Supported direct PDBQT manipulation without file format conversion.
* Used dynamic pointer vector to cache ligands.
* Used flyweight pattern to cache fragments.
* Supported dumping statistics and traceability of created ligands.
* Allowed users to specify the ranges of several chemical properties such as molecular weight.
* Allowed users to specify the number of failures of GA operations as a stopping criterion.
* Used docked atom coordinates to construct child ligands of the next generation.
* Parallelized mutation and crossover operations.
* Provided precompiled executables for both 32-bit and 64-bit Linux and Windows.


Contact Author
--------------

[Jacky Lee]


Logo
----

![igrow logo](https://github.com/HongjianLi/igrow/raw/master/logo.png)

Red grape is chosen as the logo for igrow because it is one of the author's favorite fruit. The logo image is collected from [Open Clip Art].


[drug design]: http://en.wikipedia.org/wiki/Drug_design
[AutoGrow]: http://autogrow.ucsd.edu
[AutoDock Vina]: http://vina.scripps.edu
[idock]: https://github.com/HongjianLi/idock
[Apache License 2.0]: http://www.apache.org/licenses/LICENSE-2.0.html
[C++11]: http://en.wikipedia.org/wiki/C++11
[Boost C++ Libraries]: http://www.boost.org
[doxygen]: http://www.doxygen.org
[Jacky Lee]: http://www.cse.cuhk.edu.hk/~hjli
[Open Clip Art]: http://www.openclipart.org

