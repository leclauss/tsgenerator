///\file main.cpp
///
///\brief Here starts the TSGenerator.
///
///This is the main.cpp file of the TSGenerator. The program uses a
///configuration file for the time series configuration.

#include <tsgenerator.hpp>
#include <tsgtypes.hpp>
#include <outputgenerator.hpp>
#include <global.hpp>
#include <iostream>
#include <vector>
#include <string>


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
int main(int argc, char *argv[]) {

  {
    tsg::par argTokens;
    tsg::par payload;

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
    int length = tsg::defaultLength;
    int window = tsg::defaultWindow;
    double delta = tsg::defaultDelta;
    double noise = tsg::defaultNoise;
    tsg::word type(tsg::defaultType);
    int size = tsg::defaultMotifSize;
    double height = tsg::defaultHeight;
    double step = tsg::defaultStep;
    int times = tsg::defaultTimes;
    double maxi = tsg::defaultMaxi;
    int smaller = tsg::defaultSmaller;
    tsg::word method(tsg::defaultMethod);
    tsg::word gen(tsg::defaultGen);
    tsg::rseqs motif;

    try {

      if (checkArg(argTokens, "-g", payload) || checkArg(argTokens,
            "--generator",
            payload)) {

        if (payload.empty()) {

          std::cerr << "ERROR: Generator is missing an argument." << std::endl;
          exit(EXIT_FAILURE);
        }

        gen = payload[0];
      }

      if (checkArg(argTokens, "-l", payload) || checkArg(argTokens, "--length",
            payload)) {

        if (payload.empty()) {

          std::cerr << "ERROR: Length is missing an argument." << std::endl;
          exit(EXIT_FAILURE);
        }

        if (payload[0].find_first_not_of("0123456789") != tsg::word::npos) {

          std::cerr << "ERROR: " << (payload[0]) << " is not a valid number!"
            << std::endl;
          exit(EXIT_FAILURE);
        }

        length = std::stoi(payload[0]);
      }

      if (checkArg(argTokens, "-w", payload) || checkArg(argTokens,
            "--window", payload)) {

        if (payload.empty()) {

          std::cerr << "ERROR: Window size is missing an argument." <<
            std::endl;
          exit(EXIT_FAILURE);
        }

        if (payload[0].find_first_not_of("0123456789") != tsg::word::npos) {

          std::cerr << "ERROR: " << (payload[0]) << " is not a valid number!"
            << std::endl;
          exit(EXIT_FAILURE);
        }

        window = std::stoi(payload[0]);
      }

      if (checkArg(argTokens, "-d", payload) || checkArg(argTokens, "--delta",
            payload)) {

        if (payload.empty()) {

          std::cerr << "ERROR: Delta is missing an argument." << std::endl;
          exit(EXIT_FAILURE);
        }

        if (!check_if_float(payload[0])) {

          std::cerr << "ERROR: Delta is not a float." << std::endl;
          exit(EXIT_FAILURE);
        }

        delta = std::stod(payload[0]);
      }

      if (checkArg(argTokens, "-no", payload) || checkArg(argTokens,
            "--noise", payload)) {

        if (payload.empty()) {

          std::cerr << "ERROR: Noise option is missing an argument." <<
            std::endl;
          exit(EXIT_FAILURE);
        }

        if (!check_if_float(payload[0])) {

          std::cerr << "ERROR: Noise argument is not a float." << std::endl;
          exit(EXIT_FAILURE);
        }

        noise = std::stod(payload[0]);
      }

      if (checkArg(argTokens, "-ty", payload) || checkArg(argTokens,
            "--type", payload)) {

        if (payload.empty()) {

          std::cerr << "ERROR: Type is missing an argument." << std::endl;
          exit(EXIT_FAILURE);
        }

        type = payload[0];
      }

      if (checkArg(argTokens, "-si", payload) || checkArg(argTokens,
            "--size", payload)) {

        if (payload.empty()) {

          std::cerr << "ERROR: Size is missing an argument." << std::endl;
          exit(EXIT_FAILURE);
        }

        if (payload[0].find_first_not_of("0123456789") != tsg::word::npos) {

          std::cerr << "ERROR: " << (payload[0]) << " is not a valid number!"
            << std::endl;
          exit(EXIT_FAILURE);
        }

        size = std::stoi(payload[0]);
      }

      if (checkArg(argTokens, "-he", payload) || checkArg(argTokens,
            "--height", payload)) {

        if (payload.empty()) {

          std::cerr << "ERROR: Height is missing an argument." << std::endl;
          exit(EXIT_FAILURE);
        }

        if (!check_if_float(payload[0])) {

          std::cerr << "ERROR: " << payload[0] << " is not a valid float!" <<
            std::endl;
          throw(EXIT_FAILURE);
        }

        height = std::stod(payload[0]);
      }

      if (checkArg(argTokens, "-st", payload) || checkArg(argTokens,
            "--step", payload)) {

        if (payload.empty()) {

          std::cerr << "ERROR: Step is missing an argument." << std::endl;
          exit(EXIT_FAILURE);
        }

        if (!check_if_float(payload[0])) {

          std::cerr << "ERROR: " << payload[0] << " is not a valid float!" <<
            std::endl;
          throw(EXIT_FAILURE);
        }

        step = std::stod(payload[0]);
      }

      if (checkArg(argTokens, "-ti", payload) || checkArg(argTokens,
            "--times", payload)) {

        if (payload.empty()) {

          std::cerr << "ERROR: Times is missing an argument." << std::endl;
          exit(EXIT_FAILURE);
        }

        if (payload[0].find_first_not_of("0123456789") != tsg::word::npos) {

          std::cerr << "ERROR: " << (payload[0]) << " is not a valid number!"
            << std::endl;
          exit(EXIT_FAILURE);
        }

        times = std::stoi(payload[0]);
      }

      if (checkArg(argTokens, "-me", payload) || checkArg(argTokens,
            "--Method", payload)) {

        if (payload.empty()) {

          std::cerr << "ERROR: Type is missing an argument." << std::endl;
          exit(EXIT_FAILURE);
        }

        method = payload[0];
      }

      if (checkArg(argTokens, "-ma", payload) || checkArg(argTokens,
            "--maximum", payload)) {

        if (payload.empty()) {

          std::cerr << "ERROR: Maximum is missing an argument." << std::endl;
          exit(EXIT_FAILURE);
        }

        if (!check_if_float(payload[0])) {

          std::cerr << "ERROR: " << payload[0] << " is not a valid float!" <<
            std::endl;
          throw(EXIT_FAILURE);
        }

        maxi = std::stod(payload[0]);
      }

      if (checkArg(argTokens, "-sm", payload) || checkArg(argTokens,
            "--smaller", payload)) {

        if (payload.empty()) {

          std::cerr << "ERROR: Number of smaller motifs is missing an " <<
            "argument." << std::endl;
          exit(EXIT_FAILURE);
        }

        if (payload[0].find_first_not_of("0123456789") != tsg::word::npos) {

          std::cerr << "ERROR: " << (payload[0]) << " is not a valid number!"
            << std::endl;
          exit(EXIT_FAILURE);
        }

        smaller = std::stoi(payload[0]);
      }

      //generate the time series
      tsg::rseq timeSeries;
      tsg::rseq dVector;
      tsg::iseqs motifPositions;
      bool success = false;

      try {

        tsg::TSGenerator tSGenerator(length, window, delta, noise, type, size,
            height, step, times, method, maxi, gen, smaller);
        tSGenerator.run(timeSeries, motif, dVector, motifPositions);

        success = true;
      }
      catch (...) {

        std::cerr << "ERROR: Cannot generate time series!" << std::endl;
        exit(EXIT_FAILURE);
      }

      if (success) {

        //output stuff
        OutputGenerator outputFile;

        if (checkArg(argTokens, "-o", payload) || checkArg(argTokens, "--out",
              payload)) {

          if (payload.empty()) {

            std::cerr << "ERROR: The output file name is unset." << std::endl;
            exit(EXIT_FAILURE);
          }
          else
            outputFile.setFileName(payload[0]);
        }

        outputFile.open();

        if (checkArg(argTokens, "-ho", payload) || checkArg(argTokens,
              "--horizontalOutput", payload)) {

          outputFile.printTimeSeriesHorizontal(timeSeries);
        }
        else {

          outputFile.printTimeSeriesVertical(timeSeries);
        }

        if (gen == "set motif") {

          outputFile.printMetaLine((tsg::iseq){ window }, "window");
          outputFile.printMetaLine((tsg::iseq){ motifPositions[0][0] },
              "set motif");
          outputFile.printMetaLine(motifPositions[0], "matchings");
          outputFile.printMetaLine((tsg::rseq){ dVector[0] }, "range");
          outputFile.printMetaLine(motifPositions[1], "pair motif");
          outputFile.printMetaLine((tsg::rseq){ dVector[1] }, "pair "
              "motif distance");
        }
        else if (gen == "latent motif") {

          outputFile.printMetaLine((tsg::iseq){ window }, "window");
          outputFile.printMetaLine(motif[0], "latent motif");
          outputFile.printMetaLine(motifPositions[0], "matchings");
          outputFile.printMetaLine((tsg::rseq){ dVector[0] }, "range");
          outputFile.printMetaLine(motifPositions[1], "pair motif");
          outputFile.printMetaLine((tsg::rseq){ dVector[1] }, "pair "
              "motif distance");
        }
        else {

          outputFile.printMetaLine((tsg::iseq){ window }, "window");
          outputFile.printMetaLine(motifPositions[0], "pair motif");
          outputFile.printMetaLine(dVector, "distance");
        }

        outputFile.close();
      }
    }
    catch (...) {

      std::cerr << "ERROR: Something unexpected happened!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  exit(EXIT_SUCCESS);
}
