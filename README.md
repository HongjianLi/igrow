igrow
=====

igrow is an automatic structure-based [drug design] tool for computer-aided drug discovery. It was inspired by [AutoGrow], and is hosted on GitHub at https://github.com/HongjianLi/igrow under [Apache License 2.0].


Features
--------

* igrow uses [idock] as the backend docking engine.
* igrow implements branch exchange as its crossover operator.
* igrow reads and writes ligands in pdbqt format without file format conversion.
* igrow uses an io service pool to reuse threads and maintain a high CPU utilization throughout the entire synthsizing procedure.
* igrow utilizes a dynamic pointer vector to cache and sort generated ligands.
* igrow traces the sources of generated ligands and outputs the statistics in csv format.


Supported operating systems and compilers
-----------------------------------------

* Arch Linux x86_64 and clang 3.7.0
* Mac OS X x86_64 and clang 3.7.0
* Windows 8.1 x64 and Visual Studio 2015 Update 1

Statically compiled 64-bit executables can be found in the `bin` directory.


Compilation from source code
----------------------------

igrow depends on [Boost C++ Libraries]. Boost 1.59.0 was tested. The must-be-built libraries required by igrow are `System`, `Filesystem` and `Program Options`. An unofficial header-only library, Boost.Process, is also required by igrow. The file `process.zip` must be extracted to the Boost distribution tree in order to pass compilation.

### Compilation on Linux

The Makefile uses clang as the default compiler. To compile, simply run

    make

One may modify the Makefile to use a different compiler or different compilation options.

The generated objects will be placed in the `obj` folder, and the generated executable will be placed in the `bin` folder.

### Compilation on Windows

Visual Studio 2015 solution and project files are provided. To compile, simply run

    msbuild /t:Build /p:Configuration=Release

Or one may open `igrow.sln` in Visual Studio 2015 and do a full rebuild.

The generated objects will be placed in the `obj` folder, and the generated executable will be placed in the `bin` folder.


Usage
-----

First add igrow to your PATH environment variable.

To display a full list of available options, simply run the program without arguments

    igrow

The `examples` folder contains several use cases. For example, to grow Staurosporine and dock the generated ligands against CDK2 of PDB ID 1AQ1,

    cd examples/1AQ1

One can supply the options from command line arguments

    igrow --idock_example ../../../idock/examples/1AQ1

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

* Used idock as the backend docking engine.
* Supported direct PDBQT manipulation without file format conversion.
* Used a dynamic pointer vector to cache ligands.
* Supported outputing statistics and traceability of created ligands.
* Used docked atom coordinates to construct child ligands of the next generation.
* Parallelized crossover operations.
* Provided precompiled executables for 64-bit Linux, Mac and Windows.


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
[Boost C++ Libraries]: http://www.boost.org
[doxygen]: http://www.doxygen.org
[Jacky Lee]: http://www.cse.cuhk.edu.hk/~hjli
