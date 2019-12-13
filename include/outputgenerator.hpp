///\file outputgenerator.hpp
///
///\brief File contains the output generator class declaration.
///
///This is the header file of the output generator. The output generator writes
///the time series and the motif positions into output files.
///

#ifndef OUTPUTGENERATOR_HPP
#define OUTPUTGENERATOR_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <list>
#include <cmath>
#include <cfloat>
#include <cerrno>
#include <experimental/filesystem>

using namespace std;
using namespace std::experimental::filesystem;


// --------------------------------------------------------------------------
///\brief This class represents a output generator.
///
///The OutputGenerator class offers a horizontal and vertical data print
///functions. Therefore, The time series and the motif positions list will be
///printed either horizontal or vertical into the output files.
// --------------------------------------------------------------------------
class OutputGenerator {

protected:

  // --------------------------------------------------------------------------
  ///\brief Stores the time series output file pointer.
  ///
  ///This variable stores the pointer pointing at the time series output file.
  // --------------------------------------------------------------------------
  ofstream outputFile;

  // --------------------------------------------------------------------------
  ///\brief Stores the motif sets output file pointer.
  ///
  ///This variable stores the pointer pointing at the output file containing
  ///the motif sets distances and positions.
  // --------------------------------------------------------------------------
  ofstream outputFileMotifSets;

  // --------------------------------------------------------------------------
  ///\brief Stores the GNUPlot script output file pointer.
  ///
  ///This variable stores the pointer pointing at the output file containing
  ///the GNUPlot script.
  // --------------------------------------------------------------------------
  ofstream outputFileGNUPlotScript;

  // --------------------------------------------------------------------------
  ///\brief Stores the name of the motif sets output file.
  ///
  ///This variable stores the name of the file containing the motif sets
  ///distances and positions.
  // --------------------------------------------------------------------------
  string outputFileMotifSetsName;

  // --------------------------------------------------------------------------
  ///\brief Stores the name of the GNUPlot script output file.
  ///
  ///This variable stores the name of the file containing the GNUPlot script.
  // --------------------------------------------------------------------------
  string outputFileGNUPlotScriptName;

  // --------------------------------------------------------------------------
  ///\brief Stores the prefix of the time series file.
  ///
  ///This variable stores the prefix of the time series output file.
  // --------------------------------------------------------------------------
  string outputFileName;

  // --------------------------------------------------------------------------
  ///\brief Stores the time series file name.
  ///
  ///This variable stores the time series output file name.
  // --------------------------------------------------------------------------
  string basicOutputFileName;

  // --------------------------------------------------------------------------
  ///\brief Stores the name of the time series.
  ///
  ///This variable stores the name of the time series used for labeling in
  ///gnuplot.
  // --------------------------------------------------------------------------
  string timeSeriesName = "Time Series";

  // --------------------------------------------------------------------------
  ///\brief Stores the number of the output file.
  ///
  ///This variable stores the nubmer of the current output file.
  // --------------------------------------------------------------------------
  int outputFolderNumber = 0;

  // --------------------------------------------------------------------------
  ///\brief Opens files for the output of the time series and the motif set
  ///position list.
  ///
  ///This function opens the files for the time series and the motif set
  ///position list.
  // --------------------------------------------------------------------------
  void openFile();

public:

  // --------------------------------------------------------------------------
  ///\brief Initializes the output generator.
  ///
  ///The constructor creates a new output directory path and sets up the file
  ///name with default file name "time_series". It chooses another directory
  ///path, if a directory with the same name already exists.
  // --------------------------------------------------------------------------
  OutputGenerator();

  // --------------------------------------------------------------------------
  ///\brief Initializes the output generator, i.e. creates output directory.
  ///
  ///\param [in] outputFileName_in Hands over the output file name.
  ///
  ///The constructor creates a new output directory path and sets up the file
  ///name with default file name "time_series". It chooses another directory
  ///path, if a directory with the same name already exists.
  // --------------------------------------------------------------------------
  OutputGenerator(const string &outputFileName_in);

  // --------------------------------------------------------------------------
  ///\brief Destroys the file pointers.
  ///
  ///The destructor destroys the time series file pointer and the motif set
  ///positions file pointer.
  // --------------------------------------------------------------------------
      ~OutputGenerator();

  // --------------------------------------------------------------------------
  ///\brief Writes the data horizontal into the output file.
  ///
  ///\param [in] &timeSeries_in Hands over the time series that output
  ///generator writes into output file.
  ///\param [in] &d_in Hands over the motif set distances.
  ///\param [in] &motifSetPositions_in Hands over the motif set positions that
  ///output generator writes into the output file.
  ///
  ///This function writes the time series and the motif set distances as well
  ///as positions horizontal into the output file. For that, it writes them
  ///line by line into the output files.
  // --------------------------------------------------------------------------
  void printTimeSeriesHorizontal(const vector<double> &timeSeries_in, const
      vector<double> &d_in, const vector<int> &windowSizes_in, const
      vector<vector<int>> &motifSetPositions_in);

  // --------------------------------------------------------------------------
  ///\brief Writes the data vertical into the output files.
  ///
  ///\param [in] &timeSeries_in Hands over the time series that output
  ///generator writes into output file.
  ///\param [in] &d_in Hands over the motif set distances.
  ///\param [in] &motifSetPositions_in Hands over the motif set postions that
  ///output generator writes into the output file.
  ///
  ///This function writes the time series and the motif set distances as well
  ///as positions vertical into the output file. For that, it writes them as
  ///one column into the output files.
  // --------------------------------------------------------------------------
  void printTimeSeriesVertical(const vector<double> &timeSeries_in, const
      vector<double> &d_in, const vector<int> &windowSizes_in, const
      vector<vector<int>> &motifSetPositions_in);

  // --------------------------------------------------------------------------
  ///\brief This is the name setter function for the time series label in
  ///gnuplot.
  ///
  ///\param [in] &timeSeriesName_in Hands over the name of the time series.
  ///
  ///This function sets the name of the time series label used in gnuplot.
  // --------------------------------------------------------------------------
  void setTimeSeriesName(const string &timeSeriesName_in);
};

#endif
