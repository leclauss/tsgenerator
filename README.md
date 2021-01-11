# tsgenerator

Generator of synthetic time series specialized for evaluation purposes of motif discovery algorithms.

## Introduction

This generator of synthetic time series is the result of a student research project. The idea was born, since the range of implementations generating time series specialized for evaluation of motif discovery algorithms was non-existing. This project is intended as a cross platform implementation. Therefore, the tsgenerator was tested on Ubuntu 18.04.3 LTS, Windows 10 and macOS 10.14.3 Mojave. Installation guides for each platform follow in the sections **Prerequisites** and **Installion**.

### Prerequisites

This project assumes that [Graphviz](https://www.graphviz.org/ "Graphviz website"), [Doxygen](http://doxygen.org/ "Doxygen website") and [CMake](https://cmake.org/ "CMake website") are preinstalled on the operating system. We offer a quick installation guide for everyone missing that packages in this section. 

**Linux:**

*For creating the tsgenerator-dev library*

1. Open a terminal window. Download and install the essential building tools, if they are missing.
 ```bash
 sudo apt-get install build-essential
 ```

2. Download and install [GCC](https://gcc.gnu.org/ "GNU Compiler Collection") version 8 if your gcc version is smaller.
 ```bash
 sudo apt-get install gcc-8 g++-8
 ```

3. Download and install [Git](https://git-scm.com/ "open source Website maintained by members of the Git community").
 ```bash
 sudo apt-get install git
 ```

4. Download and install [CMake](https://cmake.org/ "CMake website").
 ```bash
 sudo apt-get install cmake
 ```

5. Download and install [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/ "pkg-config website").
 ```bash
 sudo apt-get install pkg-config
 ```

*[Optional] Tools for documentation generation*

1. Download and install [Graphviz](https://www.graphviz.org/ "Graphviz website").
 ```bash
 sudo apt-get install graphviz
 ```

2. Download and install [Doxygen](http://doxygen.org/ "Doxygen website").
 ```bash
 sudo apt-get install doxygen
 ```

*[Optional] Tool for plotting with gnuplot*

1. Download and install [Gnuplot](http://gnuplot.info/ "Gnuplot website").
 ```bash
 sudo apt-get install gnuplot
 ```

*[Optional] For the GUI version of the tsgenerator (inc. plotting)*

1. Download and install [Qt5](https://www.qt.io/ "Qt website").
 ```bash
 sudo apt-get install qt5-default libqt5charts5-dev
 ```

**Windows:**

*For creating the tsgenerator-dev library*

1. Download and install [Visual Studio-IDE](https://www.visualstudio.com/
   "Microsoft Visual Studio website"). Make sure you download an IDE like Visual Studio Community, Visual Studio Professional or Visual Studio Enterprise. Visual Studio Code is an editor and not an IDE. Choose the option **Desktop development with C++** option during the installation.

2. Download and install [Git](https://git-scm.com/ "open source Website maintained by members of the Git community"). Choose the option **Use Git from Git Bash only** as well as the option **MinTTY** during installation.

3. Download and install [CMake](https://cmake.org/ "CMake website"). Choose the option **Add to system path** during the installation.

*[Optional] Tools for documentation generation*

1. Download and install [Graphviz](https://www.graphviz.org/ "Graphviz website"). Make sure you add Graphviz to the system path.

2. Download and install [Doxygen](http://doxygen.org/ "Doxygen website"). Choose the option **Add to system path** during the installation.

*[Optional] Tool for plotting with gnuplot*

1. Download and install [Gnuplot](http://gnuplot.info/ "Gnuplot website"). Make sure to hit **Add to system path** during the installation.

*[Optional] For the GUI version of the tsgenerator (inc. plotting)*

1. Download and install [Qt5](https://www.qt.io/ "Qt website"). Make sure to install the library **5.x.x** and **the QtCharts** component.

**Mac:**

*For creating the tsgenerator-dev library*

1. Download and install [Homebrew](https://brew.sh/index_de/ "Homebrew website"). Open a new terminal and run the following commands. If you downloaded a newer version of gcc change the number 8 to the appropriate version.
 ```bash
 /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
 brew install llvm gcc
 ln -s /usr/local/opt/gcc/bin/gcc-9 /usr/local/bin/gcc
 ln -s /usr/local/opt/gcc/bin/c++-9 /usr/local/bin/c++
 ln -s /usr/local/opt/llvm/bin/clang /usr/local/bin/clang
 ```

2. Download and install [Git](https://git-scm.com/ "open source website maintained by members of the Git community") and [GNU Make](https://www.gnu.org/software/make/ "GNU Make website"), if Git or GNU Make is appropriate missing.
 ```bash
 brew install git
 brew install make
 ```

3. Download and install [CMake](https://cmake.org/ "CMake website").
 ```bash
 brew install cmake
 ```

4. Download and install [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/ "pkg-config website").
 ```bash
 brew install pkg-config
 ```

*[Optional] Tools for documentation generation*

1. Download and install [Graphviz](https://www.graphviz.org/ "Graphviz website").
 ```bash
 brew install graphviz
 ```

2. Download and install [Doxygen](http://doxygen.org/ "Doxygen website").
 ```bash
 brew install doxygen
 ```

*[Optional] Tool for plotting with gnuplot*

1. Download and install [Gnuplot](http://gnuplot.info/ "Gnuplot website").
 ```bash
 brew install gnuplot
 ```

*[Optional] For the GUI version of the tsgenerator (inc. plotting)*

1. Download and install [Qt5](https://www.qt.io/ "Qt website").
 ```bash
 brew install qt5
 ```

### Installation Guide

The installation guide is divided into instructions for Linux, Windows and Mac.

**Linux:**

1. Open a terminal window. Clone the [tsgenerator](https://gitlab.com/r.moczalla/TSGenerator.git "tsgenerator project website") project directory and make the created folder the current folder.
 ```bash
 git clone "https://gitlab.com/r.moczalla/TSGenerator.git"
 cd "TSGenerator"
 ```

*Build the tsgenerator-dev library*

1. Create a build directory and make the created folder the current folder.
 ```bash
 cd tsgenerator-dev
 mkdir "build"
 cd "build"
 ```

2. Create the tsgenerator project with CMake and build the project by executing
   [GNU Make](https://www.gnu.org/software/make/ "GNU Make website"). If your default [GCC](https://gcc.gnu.org/ "GNU Compiler Collection") is at least version 8 just build the project with
 ```bash
 cmake ..
 make
 ```
 otherwise specify your [GCC](https://gcc.gnu.org/ "GNU Compiler Collection") explicitlya.
 ```bash
 CXX=g++-8 CC=gcc-8 cmake ..
 make
 ```

3. Optionally: Create debian install package.
 ```
 make package
 ```
 The package can now be installed with Ubuntus package manager.
 ```
 sudo apt install ./tsgenerator-dev-3.0.0.deb
 ```
 After you installed the package make sure to run
 ```
 sudo ldconfig
 ```

4. Optionally: Create the documentation by executing the **doc** target with GNU Make. The documentation is located in the subfolder **build/doc/**.
 ```bash
 make doc
 ```
 We offer a man page for the tsgenerator in the folder **build/doc/man/**. One may open the **TSGenerator.3** man file with the **man** command as follows.
 ```bash
 man doc/man/tsgenerator.3
 ```

To generate the tsgenerator binary repeat the steps in the directory **tsgenerator** and to generate the gui tool reapet the steps in the directory **tsgenerator-gui**.

**Windows:**

1. Open a terminal window. Clone the [tsgenerator](https://gitlab.com/r.moczalla/TSGenerator.git "tsgenerator project website") project directory and make the created folder the current folder.
 ```bash
 git clone "https://gitlab.com/r.moczalla/TSGenerator.git"
 cd "TSGenerator"
 ```

2. Create a build directory, make the created folder the current folder and build the visual studio project of tsgenerator project.
 ```bash
 mkdir "build"
 cd "build"
 cmake ..
 ```
 If you installed [Qt5](https://www.qt.io/ "Qt website") in another location than ***C:/Qt*** than specify the location ***LOC*** exclusive.
 ```bash
 cmake -DQt5_PATH="LOC" ..
 ```

3. Build the release version of the project.
 ```bash
 cmake --build . --config Release
 ```

4. Optionally: Create the documentation by building the target **doc** in the Visual Studio IDE. The documentation is located in the subfolder **build/doc/**.

To generate the tsgenerator binary repeat the steps in the directory **tsgenerator** and to generate the gui tool reapet the steps in the directory **tsgenerator-gui**.

**Mac:**

1. Open a terminal window. Clone the [tsgenerator](https://gitlab.com/r.moczalla/TSGenerator.git "tsgenerator project website") project directory and make the created folder the current folder.
 ```bash
 git clone "https://gitlab.com/r.moczalla/TSGenerator.git"
 cd "TSGenerator"
 ```

2. Create a build directory and make the created folder the current folder.
 ```bash
 mkdir "build"
 cd "build"
 ```

3. Create the tsgenerator project with CMake and build the project by executing make.
 ```bash
 cmake ..
 make
 ```

4. Optionally: Create the documentation by executing the **doc** target with make. The documentation is located in the subfolder **build/doc/**.
 ```bash
 make doc
 ```
 We offer a man page for the tsgenerator in the folder **build/doc/man/**. One may open the **TSGenerator.3** man file with the **man** command as follows.
 ```bash
 man doc/man/TSGenerator.3
 ```

To generate the tsgenerator binary repeat the steps in the directory **tsgenerator** and to generate the gui tool reapet the steps in the directory **tsgenerator-gui**.

## Running Tests

**Linux:**

Run the following command in the build directory to run the tests.
```bash
make test
```

**Windows:**

1. Build with Visual Studio the target tsgeneratorTest.

2. Open the git bash, navigate into the build directory and run the following command.
 ```bash
 ctest -C Debug
 ```

**Mac:**

Run the following command in the build directory to run the tests.
```bash
make test
```

### Quick Usage

1. One may generate a synthetic time series with a length of 10000, time series top latent motif attributes type sine, 3 top latent motif matching subseqeunces, window size 30, top latent motif height 10.0 and a randomization factor of 0.01 by running the tsgenerator with the following command. On Linux and Mac in the **build/** directory and on Windows in the **build/src/** directory.
 ```bash
 ./tsgenerator -g "latent motif" -ty sine -me "boundedNormalRandomWalk" -l 10000 -w 30 -no 1.0 -si 3 -h 10.0
 ```

2. tsgenerator command has the following options.
 * **-g STRING, --generator STRING** sets the method for injecting sequences into the time series matching the synthetic motif. The available methods are **pair motif**, **set motif** and **latent motif**.
 * **-ty STRING, --type STRING** sets the motif type, the shape of the injected motif and the inserted sequences. The available types are **box**, **triangle**, **semicircle**, **trapezoid**, **positiveflank**, **negativeflank**, **sine** and **cosine**.
 * **-me STRING, --method STRING** sets the method for generating the base time series. The available methods are **simpleRandomWalk**, **realRandomWalk**, **normalRandomWalk**, **linearRandomWalk**, **boundedSimpleRandomWalk**, **boundedRealRandomWalk**, **boundedNormalRandomWalk**, **boundedLinearRandomWalk**, **uniformRandom**, **normalRandom**, **piecewiseLinearRandom** and **splineRepeated**.
 * **-l INTEGER, --length INTEGER** sets the length of the synthetic time series.
 * **-w INTEGER, --window INTEGER** sets the window size of the synthetic time series motif sets subsequences.
 * **-si INTEGER, --size INTEGER** sets the size of the motif, the number of inserted sequences non-self matched by the motif.
 * **-no FLOAT, --noise FLOAT** sets the noise value. A random value in the range from -FLOAT to FLOAT is added to the base time series and the inserted sequences.
 * **-d FLOAT, --delta FLOAT** sets the maximum absolute difference between two consecutive values in the time series.
 * **-he FLOAT, --height FLOAT** sets the maximum absolute difference between two values of the base motif.
 * **-st FLOAT, --step FLOAT** sets the maximum step size in x direction from two consecutive values when creating a linear approximated or splined base time series.
 * **-ti INTEGER, --times INTEGER** sets the number of values computed to generate a repeating pattern when generating a linear approsimated or splined base time series.
 * **-ma FLOAT, --maxi FLOAT** sets the maximum absolute value in the base times series.
 * **-o STRING, --out STRING** sets the base name of the output files.
 * **-ho, --horizontalOutput** prints the time series values horizontal in the output file divided by a delimiter.
 * **-h, --help** prints the help text.
 * **-v, --version** prints the version information.

### Output

The tsgenerator creates a new output folder called output and four files. A **time_series.csv** file, a **time_series_meta.csv** file and a **time_series_plot.plt** file.
1. **time_series.csv**
 Contains the synthetic time series values.

2. **time_series_meta.csv**
 Contains the synthetic time series motif meta data. Each data set first value is the similarity of the synthetic time series top pair motif or the range of the synthetic time series top latent motif. The second value is the window size of the synthetic time series top pair motif or the synthetic time series top latent motif followed by all the positions of the subsequences in the synthetic time series top pair motif or the synthetic time series top latent motif matching subseqeunces positions in the synthetic time series.

3. **time_series_plot.plt**
 Contains a [Gnuplot](http://www.gnuplot.info/ "Gnuplot website") script file script file that plots the **time_series.csv** file.

## Built With

* [Graphviz](https://www.graphviz.org/ "Graphviz website") is an open source, cross platform graph visualization software.
* [Doxygen](http://doxygen.org/ "Doxygen website") is a cross platform documentation tool.
* [CMake](https://cmake.org/ "CMake website") is an open source, cross platform family of tools designed to build, test and package software.

## References

* [GNU C++ Library](https://gcc.gnu.org/onlinedocs/libstdc++/ "GNU C++ Library website") is a collection of functions and classes written in the core C++ language and they are part of the C++ ISO Standard.
* [Gnuplot](http://gnuplot.info/ "Gnuplot website") is an cross platform and portable command-line driven graphing utility.
* [Gnuplot i](https://gitlab.control.lth.se/letter2martin/mod_dmp_example/blob/8b719a0bc54b801b79b8f8cfd32dfbc662d930b9/gnuplot-cpp/ "Gnuplot i source website") is a cross platform C++ interface to gnuplot.

## Versioning

We version the project with [GitLab](https://gitlab.com/ "GitLab website").

## Author

* **Rafael Moczalla** - *Computer Scientist & Software Developer*

## License

MIT License

Copyright (c) 2018 Rafael Moczalla

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Acknowledgments

* To Dr. rer. nat. Patrick Schäfer for revealing discussions.
