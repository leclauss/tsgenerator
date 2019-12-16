# TSGenerator

Generator of synthetic time series specialized for evaluation purposes of motif discovery algorithms.

## Introduction

This generator of synthetic time series is the result of a student research project. The idea was born, since the range of implementations generating time series specialized for evaluation of motif discovery algorithms was non-existing. This project is intended as a cross platform implementation. Therefore, the TSGenerator was tested on Ubuntu 18.04.3 LTS, Windows 10 and macOS 10.14.3 Mojave. Installation guides for each platform follow in the sections **Prerequisites** and **Installion**.

### Prerequisites

This project assumes that [Graphviz](https://www.graphviz.org/ "Graphviz website"), [Doxygen](http://doxygen.org/ "Doxygen website") and [CMake](https://cmake.org/ "CMake website") are preinstalled on the operating system. We offer a quick installation guide for everyone missing that packages in this section. 

**Linux:**

1. Open a terminal window. Download and install the essential building tools, if they are missing.
 ```bash
 sudo apt-get install build-essential
 ```

2. Download and install [Graphviz](https://www.graphviz.org/ "Graphviz website").
 ```bash
 sudo apt-get install graphviz
 ```

3. Download and install [Doxygen](http://doxygen.org/ "Doxygen website").
 ```bash
 sudo apt-get install doxygen
 ```

4. Download and install [CMake](https://cmake.org/ "CMake website").
 ```bash
 sudo apt-get install cmake
 ```

5. Download and install [Git](https://git-scm.com/ "open source Website maintained by members of the Git community").
 ```bash
 sudo apt-get install git
 ```

6. Download and install [Gnuplot](http://gnuplot.info/ "Gnuplot website").
 ```bash
 sudo apt-get install gnuplot
 ```

7. Download and install [GCC](https://gcc.gnu.org/ "GNU Compiler Collection") version 8 if your gcc version is smaller.
 ```bash
 sudo apt-get install gcc-8 g++-8
 ```

**Windows:**

1. Download and install [Graphviz](https://www.graphviz.org/ "Graphviz website"). Make sure you add Graphviz to the system path.

2. Download and install [Doxygen](http://doxygen.org/ "Doxygen website"). Choose the option **Add to system path** during the installation.

3. Download and install [CMake](https://cmake.org/ "CMake website"). Choose the option **Add to system path** during the installation.

4. Download and install [Visual Studio-IDE](https://www.visualstudio.com/
   "Microsoft Visual Studio website"). Make sure you download an IDE like Visual Studio Community, Visual Studio Professional or Visual Studio Enterprise. Visual Studio Code is an editor and not an IDE. Choose the option **Desktop development with C++** option and in the second tab the **C++ Clang Compiler for Windows**, **C++ Clang-cl for ... build tools (x64/x86)** as well as **CMake** during the installation.

5. Download and install [Git](https://git-scm.com/ "open source Website maintained by members of the Git community"). Choose the option **Use Git from Git Bash only** as well as the option **MinTTY** during installation.

6. Download and install [Gnuplot](http://gnuplot.info/ "Gnuplot website"). Make sure to hit **Add to system path** during the installation.

**Mac:**

1. Download and install [Homebrew](https://brew.sh/index_de/ "Homebrew website"). Open a new terminal and run the following commands. If you downloaded a newer version of gcc change the number 8 to the appropriate version.
 ```bash
 /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
 brew install llvm gcc
 ln -s /usr/local/opt/gcc/bin/gcc-9 /usr/local/bin/gcc
 ln -s /usr/local/opt/gcc/bin/c++-9 /usr/local/bin/c++
 ln -s /usr/local/opt/llvm/bin/clang /usr/local/bin/clang
 ```

2. Download and install [Graphviz](https://www.graphviz.org/ "Graphviz website").
 ```bash
 brew install graphviz
 ```

3. Download and install [Doxygen](http://doxygen.org/ "Doxygen website").
 ```bash
 brew install doxygen
 ```

4. Download and install [CMake](https://cmake.org/ "CMake website").
 ```bash
 brew install cmake
 ```

5. Download and install [Git](https://git-scm.com/ "open source website maintained by members of the Git community") and [GNU Make](https://www.gnu.org/software/make/ "GNU Make website"), if Git or GNU Make is appropriate missing.
 ```bash
 brew install git
 brew install make
 ```

6. Download and install [Gnuplot](http://gnuplot.info/ "Gnuplot website").
 ```bash
 brew install gnuplot
 ```

### Installation Guide

The installation guide is divided into instructions for Linux, Windows and Mac.

**Linux:**

1. Open a terminal window. Clone the [TSGenerator](https://gitlab.com/r.moczalla/TSGenerator.git "TSGenerator project website") project directory and make the created folder the current folder.
 ```bash
 git clone "https://gitlab.com/r.moczalla/TSGenerator.git"
 cd "TSGenerator"
 ```

2. Create a build directory and make the created folder the current folder.
 ```bash
 mkdir "build"
 cd "build"
 ```

3. Create the TSGenerator project with CMake and build the project by executing
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

4. TSGenerator is located in the subfolder **bin**. Make the created subfolder the current folder in the Git Bash and run the TSGenerator.
 ```bash
 cd "bin"
 ./TSGenerator [Options]
 ```

5. Optionally: Create the documentation by executing the **doc** target with GNU Make. The documentation is located in the subfolder **build/doc/**.
 ```bash
 make doc
 ```
 We offer a man page for the TSGenerator in the folder **build/doc/man/**. One may open the **TSGenerator.3** man file with the **man** command as follows.
 ```bash
 man doc/man/TSGenerator.3
 ```

**Windows:**

1. Open a terminal window. Clone the [TSGenerator](https://gitlab.com/r.moczalla/TSGenerator.git "TSGenerator project website") project directory and make the created folder the current folder.
 ```bash
 git clone "https://gitlab.com/r.moczalla/TSGenerator.git"
 cd "TSGenerator"
 ```

2. Create a build directory, make the created folder the current folder and build the visual studio project of TSGenerator project.
 ```bash
 mkdir "build"
 cd "build"
 cmake ..
 ```

3. Open the **TSGenerator.sln** file with Visual Studio and rightclick on the target TSGenerator in the projectmap explorer. Hit the option **Select as Startproject** and press the f5 key to build the project.

4. TSGenerator is located in the subfolder **bin**. Make the created subfolder the current folder in the Git Bash and run the TSGenerator.
 ```bash
 cd "bin"
 ./TSGenerator [Options]
 ```

5. Optionally: Create the documentation by building the target **doc** in the Visual Studio IDE. The documentation is located in the subfolder **build/doc/**.

**Mac:**

1. Open a terminal window. Clone the [TSGenerator](https://gitlab.com/r.moczalla/TSGenerator.git "TSGenerator project website") project directory and make the created folder the current folder.
 ```bash
 git clone "https://gitlab.com/r.moczalla/TSGenerator.git"
 cd "TSGenerator"
 ```

2. Create a build directory and make the created folder the current folder.
 ```bash
 mkdir "build"
 cd "build"
 ```

3. Create the TSGenerator project with CMake and build the project by executing make.
 ```bash
 cmake ..
 make
 ```

4. TSGenerator is located in the subfolder **bin**. Make the created subfolder the current folder in the Git Bash and run the TSGenerator.
 ```bash
 cd "bin"
 ./TSGenerator [Options]
 ```

5. Optionally: Create the documentation by executing the **doc** target with make. The documentation is located in the subfolder **build/doc/**.
 ```bash
 make doc
 ```
 We offer a man page for the TSGenerator in the folder **build/doc/man/**. One may open the **TSGenerator.3** man file with the **man** command as follows.
 ```bash
 man doc/man/TSGenerator.3
 ```

## Running Tests

**Linux:**

Run the following command in the build directory to run the tests.
```bash
make test
```

**Windows:**

1. Build with Visual Studio the target TSGeneratorTest.

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

1. One may generate a synthetic time series with a length of 10000, time series top latent motif attributes type sine, 3 top latent motif matching subseqeunces, window size 50, top latent motif height 100.0 and a randomization factor of 0.01 by running the TSGenerator with the following command. On Linux and Mac in the **build/** directory and on Windows in the **build/src/** directory.
 ```bash
 ./TSGenerator -l 10000 -w 50 -rd 0.01 -lm sine 3 100.0
 ```
 
2. TSGenerator command has the following syntax.
 ```bash
 ./TSGenerator -l INTEGER -w INTEGER -rd FLOAT [Options]
 ```
 
 The argument **-l** sets the length of the synthetic time series.
 The argument **-w** sets the window size of the synthetic time series motif sets subsequences.
 The argument **-rd FLOAT** sets the randomization factor that is used to generate the sythetic time series values as well as the synthetic time series top pair motifs or the time series top latent motif matching subsequences.
 
 The available options are -db, -pm, -lm, -o, -tsn, -ho, -r, -rd, -h and -v.
 * **-db** sets the operating mode to data base.
 * **-pm** sets the operating mode to top pair motif and the subsequences type as well as height. The default operating mode is the top pair motif mode, the default type is BOX and the default height is 200. The synthax is STRING FLOAT.
 * **-lm** sets the operating mode to top latent motif and the synthetic time series top latent motif attributes. The attributes are the type, the count and the height of the synthetic time series top latent motif. The synthax is STRING INTEGER FLOAT.
 * **-o STRING** sets the output file name.
 * **-tsn STRING** sets the synthetic time series name.
 * **-ho** sets the output to horizontal mode. Each row in the CSV files is now a data set.
 * **-r FLOAT FLOAT** sets the range of the synthetic time series values.
 * **-h** prints the help text.
 * **-v** prints the version information.

### Output

The TSGenerator creates a new output folder called output and four files. A **time_series.csv** file, a **time_series_meta.csv** file and a **time_series_plot.plt** file. 
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

Copyright (c) 2019 Rafael Moczalla

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Acknowledgments

* To Dr. rer. nat. Patrick Schäfer for revealing discussions.
