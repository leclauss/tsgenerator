///\file main.cpp
///
///\brief Here starts the TSGenerator.
///
///This is the main.cpp file of the TSGenerator. The program uses a
///configuration file for the time series configuration.

#ifdef _WIN32
#include <Windows.h>
#endif

#include <tsgenerator.hpp>
#include <tsgtypes.hpp>
#include <outputgenerator.hpp>
#include <global.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <QApplication>
#include <QWidget>


using namespace std;


// ----------------------------------------------------------------------------
///\brief Function called at program start.
///
///\param [in] argc Hands over the number of arguments to main.
///\param [in] **argv Hands over the arguments to main.
///
///This is the typical main function of C++. The main function is called at the
///program start. First all non-local objects with static storage duration are
///initialized. We use the main function to parse the argument list to set up
///and to start the TSGenerator. Also a configuration file is created, if
///a configuration file do not exist.
// ----------------------------------------------------------------------------

int main(int argc, char *argv[]) {

#ifdef _WIN32
  FreeConsole();
#endif

  {
    //default parameters
    int length = 1000;
    int window = 30;
    double delta = 1.0;
    double noise = 0.1;
    string type("box");
    int size = 3;
    double height = 10.0;
    double step = 1.0;
    int times = 3;
    double maxi = 20.0;
    string method("boundedNormalRandomWalk");

    try {

      //generate the time series
      tsg::rseq timeSeries;
      tsg::rseq dVector;
      tsg::iseq windows;
      tsg::iseqs motifPositions;
      tsg::TSGenerator tSGenerator(length, window, delta, noise, type, size,
          height, step, times, method, maxi);
      tSGenerator.run(timeSeries, dVector, windows, motifPositions);

      //output stuff
      OutputGenerator *outputFile = nullptr;

      outputFile = new OutputGenerator();

      outputFile->printTimeSeriesVertical(timeSeries, dVector, windows,
          motifPositions);

      delete outputFile;
    }
    catch (int e) {

      if (e == EXIT_FAILURE)
        exit(e);
      else {

        cerr << "ERROR: Something unexpected happened!" << endl;
        exit(EXIT_FAILURE);
      }
    }
    catch (...) {

      cerr << "ERROR: Something unexpected happened!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  QApplication app(argc, argv);

  QWidget window;

  window.resize(250, 150);
  window.setWindowTitle("tsgenerator");
  window.show();

  return app.exec();
}
