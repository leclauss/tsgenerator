///\file main.cpp
///
///\brief Here starts the TSGenerator.
///
///This is the main.cpp file of the TSGenerator. The program uses a
///configuration file for the time series configuration.

#include <tsgenerator.hpp>
#include <outputgenerator.hpp>
#include <global.hpp>
#include <iostream>
#include <vector>
#include <string>

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

  {
    vector<string> argTokens;
    vector<string> payload;

    parseArgs(argc, argv, argTokens);

    if (checkArg(argTokens, "-h", payload) || checkArg(argTokens, "--help",
          payload)) {

      print_help();
      exit(EXIT_SUCCESS);
    }

    if (checkArg(argTokens, "-v", payload) || checkArg(argTokens, "--version",
         payload)) {

      print_version();
      exit(EXIT_SUCCESS);
    }

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

      if (checkArg(argTokens, "-l", payload) || checkArg(argTokens, "--length",
            payload)) {

        if (payload.empty()) {

          cerr << "ERROR: Length is missing an argument." << endl;
          exit(EXIT_FAILURE);
        }

        if (payload[0].find_first_not_of("0123456789") != string::npos) {

          cerr << "ERROR: " << (payload[0]) << " is not a valid number!" << endl;
          exit(EXIT_FAILURE);
        }

        length = stoi(payload[0]);
      }

      if (checkArg(argTokens, "-w", payload) || checkArg(argTokens,
            "--window", payload)) {

        if (payload.empty()) {

          cerr << "ERROR: Window size is missing an argument." << endl;
          exit(EXIT_FAILURE);
        }

        if (payload[0].find_first_not_of("0123456789") != string::npos) {

          cerr << "ERROR: " << (payload[0]) << " is not a valid number!" << endl;
          exit(EXIT_FAILURE);
        }

        window = stoi(payload[0]);
      }

      if (checkArg(argTokens, "-d", payload) || checkArg(argTokens, "--delta",
            payload)) {

        if (payload.empty()) {

          cerr << "ERROR: Delta is missing an argument." << endl;
          exit(EXIT_FAILURE);
        }

        if (!check_if_float(payload[0])) {

          cerr << "ERROR: Delta is not a float." << endl;
          exit(EXIT_FAILURE);
        }

        delta = stod(payload[0]);
      }

      if (checkArg(argTokens, "-no", payload) || checkArg(argTokens,
            "--noise", payload)) {

        if (payload.empty()) {

          cerr << "ERROR: Noise option is missing an argument." << endl;
          exit(EXIT_FAILURE);
        }

        if (!check_if_float(payload[0])) {

          cerr << "ERROR: Noise argument is not a float." << endl;
          exit(EXIT_FAILURE);
        }

        noise = stod(payload[0]);
      }

      if (checkArg(argTokens, "-ty", payload) || checkArg(argTokens,
            "--type", payload)) {

        if (payload.empty()) {

          cerr << "ERROR: Type is missing an argument." << endl;
          exit(EXIT_FAILURE);
        }

        type = payload[0];
      }

      if (checkArg(argTokens, "-si", payload) || checkArg(argTokens,
            "--size", payload)) {

        if (payload.empty()) {

          cerr << "ERROR: Size is missing an argument." << endl;
          exit(EXIT_FAILURE);
        }

        if (payload[0].find_first_not_of("0123456789") != string::npos) {

          cerr << "ERROR: " << (payload[0]) << " is not a valid number!" << endl;
          exit(EXIT_FAILURE);
        }

        size = stoi(payload[0]);
      }

      if (checkArg(argTokens, "-he", payload) || checkArg(argTokens,
            "--height", payload)) {

        if (payload.empty()) {

          cerr << "ERROR: Height is missing an argument." << endl;
          exit(EXIT_FAILURE);
        }

        if (!check_if_float(payload[0])) {

          cerr << "ERROR: " << payload[0] << " is not a valid float!" << endl;
          throw(EXIT_FAILURE);
        }

        height = stoll(payload[0]);
      }

      if (checkArg(argTokens, "-st", payload) || checkArg(argTokens,
            "--step", payload)) {

        if (payload.empty()) {

          cerr << "ERROR: Step is missing an argument." << endl;
          exit(EXIT_FAILURE);
        }

        if (!check_if_float(payload[0])) {

          cerr << "ERROR: " << payload[0] << " is not a valid float!" << endl;
          throw(EXIT_FAILURE);
        }

        step = stoll(payload[0]);
      }

      if (checkArg(argTokens, "-ti", payload) || checkArg(argTokens,
            "--times", payload)) {

        if (payload.empty()) {

          cerr << "ERROR: Times is missing an argument." << endl;
          exit(EXIT_FAILURE);
        }

        if (payload[0].find_first_not_of("0123456789") != string::npos) {

          cerr << "ERROR: " << (payload[0]) << " is not a valid number!" << endl;
          exit(EXIT_FAILURE);
        }

        times = stoi(payload[0]);
      }

      if (checkArg(argTokens, "-me", payload) || checkArg(argTokens,
            "--Method", payload)) {

        if (payload.empty()) {

          cerr << "ERROR: Type is missing an argument." << endl;
          exit(EXIT_FAILURE);
        }

        method = payload[0];
      }

      if (checkArg(argTokens, "-ma", payload) || checkArg(argTokens,
            "--maximum", payload)) {

        if (payload.empty()) {

          cerr << "ERROR: Maximum is missing an argument." << endl;
          exit(EXIT_FAILURE);
        }

        if (!check_if_float(payload[0])) {

          cerr << "ERROR: " << payload[0] << " is not a valid float!" << endl;
          throw(EXIT_FAILURE);
        }

        maxi = stoll(payload[0]);
      }

      //generate the time series
      vector<double> timeSeries;
      vector<double> dVector;
      vector<int> windows;
      vector<vector<int>> motifPositions;
      TSGenerator tSGenerator(length, window, delta, noise, type, size, height,
          step, times, method, maxi);
      tSGenerator.run(timeSeries, dVector, windows, motifPositions);

      //output stuff
      OutputGenerator *outputFile = nullptr;

      if (checkArg(argTokens, "-o", payload) || checkArg(argTokens, "--out",
            payload)) {

        if (payload.empty()) {

          cerr << "ERROR: The output file name is unset." << endl;
          exit(EXIT_FAILURE);
        }
        else
          outputFile = new OutputGenerator(payload[0]);
      }
      else
        outputFile = new OutputGenerator();

      if (checkArg(argTokens, "-tsn", payload) || checkArg(argTokens,
            "--timeSeriesName", payload)) {

        if (payload.empty()) {

          cerr << "ERROR: The time series name is unset." << endl;
          exit(EXIT_FAILURE);
        }
        else
          outputFile->setTimeSeriesName(payload[0]);
      }

      if (checkArg(argTokens, "-ho", payload) || checkArg(argTokens,
            "--horizontalOutput", payload))
        outputFile->printTimeSeriesHorizontal(timeSeries, dVector, windows,
            motifPositions);
      else
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

  exit(EXIT_SUCCESS);
}
