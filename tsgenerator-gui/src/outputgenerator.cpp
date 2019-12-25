///\file outputgenerator.cpp
///
///\brief File contains the output generator class definition.
///
///This is the source file of the output generator. The output generator writes
///the time series and the motif set positions into output files.

#include <outputgenerator.hpp>

OutputGenerator::OutputGenerator() : OutputGenerator("time_series") {}

OutputGenerator::OutputGenerator(const string &outputFileName_in)
  : outputFileName("/" + outputFileName_in) {

  string outputDirName = "./time_series_";

  while (exists(outputDirName + to_string(outputFolderNumber)))
    outputFolderNumber++;

  create_directory(outputDirName + to_string(outputFolderNumber));

  basicOutputFileName = outputFileName;
  outputFileMotifSetsName = outputDirName + to_string(outputFolderNumber)
    + outputFileName + "_meta";
  outputFileGNUPlotScriptName = outputDirName + to_string(outputFolderNumber)
    + outputFileName + "_plot";
  outputFileName = outputDirName + to_string(outputFolderNumber)
    + outputFileName;

  basicOutputFileName.erase(0, 1);
}

OutputGenerator::~OutputGenerator() {

  if (outputFile.is_open())
    outputFile.close();

  if (outputFileMotifSets.is_open())
    outputFileMotifSets.close();

  if (outputFileGNUPlotScript.is_open())
    outputFileGNUPlotScript.close();
}

void OutputGenerator::openFile() {

  if (outputFile.is_open())
    outputFile.close();

  if (outputFileMotifSets.is_open())
    outputFileMotifSets.close();

  if (outputFileGNUPlotScript.is_open())
    outputFileGNUPlotScript.close();

  string tmpTimeSeriesFileName = outputFileName + "_"
    + to_string(outputFolderNumber) + ".csv";
  string tmpFileMotifSetPositionsName = outputFileMotifSetsName + "_"
    + to_string(outputFolderNumber) + ".csv";
  string tmpGNUPlotScriptName = outputFileGNUPlotScriptName + "_"
    + to_string(outputFolderNumber) + ".plt";

  if (exists(tmpTimeSeriesFileName)) {

    cerr << "ERROR: Time series file already exist!" << endl;
    throw(EXIT_FAILURE);
  }

  if (exists(tmpFileMotifSetPositionsName)) {

    cerr << "ERROR: Motif positions file already exist!" << endl;
    throw(EXIT_FAILURE);
  }

  if (exists(tmpGNUPlotScriptName)) {

    cerr << "ERROR: GNUPlot script file already exist!" << endl;
    throw(EXIT_FAILURE);
  }

  outputFile.open(tmpTimeSeriesFileName, ios::out);
  outputFileMotifSets.open(tmpFileMotifSetPositionsName, ios::out);
  outputFileGNUPlotScript.open(tmpGNUPlotScriptName, ios::out);
}

void OutputGenerator::printTimeSeriesHorizontal(const vector<double>
    &timeSeries_in, const vector<double> &d_in, const vector<int>
    &windowSizes_in, const vector<vector<int>> &motifSetPositions_in) {

  if (timeSeries_in.size() == 0) {

    cerr << "ERROR: Time series is empty." << endl;
    throw(EXIT_FAILURE);
  }

  if (d_in.size() == 0) {

    cerr << "ERROR: Distance vector is empty." << endl;
    throw(EXIT_FAILURE);
  }

  if (d_in.size() != motifSetPositions_in.size()) {

    cerr << "ERROR: Distance vector and motif set positions vector" <<
      " have different size." << endl;
    throw(EXIT_FAILURE);
  }

  if (windowSizes_in.size() == 0) {

    cerr << "ERROR: Window size vector is empty." << endl;
    throw(EXIT_FAILURE);
  }

  if (windowSizes_in.size() != d_in.size()) {

    cerr << "ERROR: Window size vector and distance vector have" <<
      " different size." << endl;
    throw(EXIT_FAILURE);
  }

  //check if motif set positions vector is empty
  int maxMSSize = 0;

  for (int motifSetItr = 0; motifSetItr < (int)d_in.size(); motifSetItr++) {

    if ((motifSetPositions_in[motifSetItr]).size() == 0) {

      cerr << "ERROR: Motif set positions vector is empty." << endl;
      throw(EXIT_FAILURE);
    }
    else if ((int)(motifSetPositions_in[motifSetItr]).size() > maxMSSize)
      maxMSSize = (motifSetPositions_in[motifSetItr]).size();
  }


  openFile();

  //print the time series
  int itr = 1;
  double valueMin = DBL_MAX, valueMax = DBL_MIN;
  double tmpOutputValue;
  string tmpOutputString;

  tmpOutputValue = timeSeries_in[0];
  tmpOutputString = to_string(tmpOutputValue);
  outputFile.write(tmpOutputString.c_str(), tmpOutputString.size());

  while (itr < (int)timeSeries_in.size()) {

    tmpOutputValue = timeSeries_in[itr];
    tmpOutputString = ", " + to_string(tmpOutputValue);
    outputFile.write(tmpOutputString.c_str(), tmpOutputString.size());

    if (tmpOutputValue > valueMax)
      valueMax = tmpOutputValue;
    if (tmpOutputValue < valueMin)
      valueMin = tmpOutputValue;

    itr++;
  }

  outputFile.close();


  tmpOutputString = "\n";

  string gnuplotCMD;
#ifdef _WIN32
  gnuplotCMD = "set terminal wxt size 1600, 300\n";
#elif __APPLE__
  gnuplotCMD = "set terminal aqua size 1600, 300\n";
#else
  gnuplotCMD = "set terminal qt size 1600, 300 persist\n";
#endif
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
  gnuplotCMD = "set xrange[" + to_string(0) + ":" + to_string(itr - 1) + "]\n";
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
  gnuplotCMD = "set yrange[" + to_string(valueMin - 0.2 * (valueMax
        - valueMin)) + ":" + to_string(valueMax + 0.2 * (valueMax - valueMin))
    + "]\n";
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
  gnuplotCMD = "set grid\n";
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());


  int maxMotifSetSize = 0;

  for (vector<int> item : motifSetPositions_in)
    if ((int)item.size() > maxMotifSetSize)
      maxMotifSetSize = item.size();

  tmpOutputString = " , \"range/similarity\"";
  outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());

  tmpOutputString = ", \"window size\"";
  outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());

  for (int vectorItr = 0; vectorItr < maxMotifSetSize; vectorItr++) {

    tmpOutputString = ", \"position " + to_string(vectorItr) + "\"";
    outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());
  }

  tmpOutputString = "\n";
  outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());

  for (int vectorItr = 0; vectorItr < (int)d_in.size(); vectorItr++) {

    if (motifSetPositions_in[vectorItr].size() == 2)
      tmpOutputString = "\"Top Pair Motif Locations\"";
    else
      tmpOutputString = "\"Top Latent Motif Locations\"";
    outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());

          ostringstream streamObj;
          streamObj << fixed << setprecision(6) << ceil(1000000.0
              * d_in[vectorItr]) / 1000000.0;

    tmpOutputString = ", " + streamObj.str();
    outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());

    tmpOutputString =  ", " + to_string(windowSizes_in[vectorItr]);
    outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());

    for (int positionItr = 0; positionItr
        < (int)(motifSetPositions_in[vectorItr]).size(); positionItr++) {

      tmpOutputString = ", "
        + to_string((motifSetPositions_in[vectorItr])[positionItr]);
      outputFileMotifSets.write(tmpOutputString.c_str(),
          tmpOutputString.size());

      if (motifSetPositions_in[vectorItr].size() != 2) {

        gnuplotCMD = "set obj rect from "
          + to_string(motifSetPositions_in[vectorItr][positionItr])
          + ", graph 0 to "
          + to_string(motifSetPositions_in[vectorItr][positionItr]
              + windowSizes_in[vectorItr] - 1)
          + ", graph 1 back fc rgb \"#c5c5c5\" "
          + "fs border rgb \"#c5c5c5\"\n";
        outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
      }
    }

    tmpOutputString = "\n";
    outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());
  }

  for (int vectorItr = 0; vectorItr < (int)d_in.size(); vectorItr++) {

    for (int positionItr = 0; positionItr
        < (int)(motifSetPositions_in[vectorItr]).size(); positionItr++) {

      if (motifSetPositions_in[vectorItr].size() == 2) {

        gnuplotCMD = "set obj rect from "
          + to_string(motifSetPositions_in[vectorItr][positionItr])
          + ", graph 0 to "
          + to_string(motifSetPositions_in[vectorItr][positionItr]
              + windowSizes_in[vectorItr] - 1)
          + ", graph 1 back fc rgb \"#656565\" "
          + "fs pattern 1 border rgb \"#656565\"\n";
        outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
      }
    }
  }


  gnuplotCMD = "plot \"" + basicOutputFileName + "_"
    + to_string(outputFolderNumber) + ".csv\" matrix title \"" + timeSeriesName
    + " " + to_string(outputFolderNumber) + "\" with lines";
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());

#ifdef _WIN32
#elif __APPLE__
#else
  gnuplotCMD = "\n\npause mouse close";
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
#endif
  outputFileGNUPlotScript.close();

  outputFileMotifSets.close();
}

void OutputGenerator::printTimeSeriesVertical(const vector<double>
    &timeSeries_in, const vector<double> &d_in, const vector<int>
    &windowSizes_in, const vector<vector<int>> &motifSetPositions_in) {

  if (timeSeries_in.size() == 0) {

    cerr << "ERROR: Time series is empty." << endl;
    throw(EXIT_FAILURE);
  }

  if (d_in.size() == 0) {

    cerr << "ERROR: Distance vector is empty." << endl;
    throw(EXIT_FAILURE);
  }

  if (d_in.size() != motifSetPositions_in.size()) {

    cerr << "ERROR: Distance vector and motif set positions vector" <<
      " have different size." << endl;
    throw(EXIT_FAILURE);
  }

  if (windowSizes_in.size() == 0) {

    cerr << "ERROR: Window size vector is empty." << endl;
    throw(EXIT_FAILURE);
  }

  if (windowSizes_in.size() != d_in.size()) {

    cerr << "ERROR: Window size vector and distance vector have" <<
      " different size." << endl;
    throw(EXIT_FAILURE);
  }

  //check if there is an empty motif set vector
  int maxMSSize = 0;

  for (int motifSetItr = 0; motifSetItr < (int)d_in.size(); motifSetItr++) {

    if ((int)(motifSetPositions_in[motifSetItr]).size() == 0) {

      cerr << "ERROR: Motif set positions vector is empty." << endl;
      throw(EXIT_FAILURE);
    }
    else if ((int)(motifSetPositions_in[motifSetItr]).size() > maxMSSize)
      maxMSSize = (int)(motifSetPositions_in[motifSetItr]).size();
  }


  openFile();

  //print the time series
  int itr = 0;
  double valueMin = DBL_MAX, valueMax = DBL_MIN;
  double tmpOutputValue;
  string tmpOutputString;

  while (itr < (int)timeSeries_in.size()) {

    tmpOutputValue = timeSeries_in[itr];
    tmpOutputString = to_string(tmpOutputValue) + "\n";
    outputFile.write(tmpOutputString.c_str(), tmpOutputString.size());

    if (tmpOutputValue > valueMax)
      valueMax = tmpOutputValue;
    if (tmpOutputValue < valueMin)
      valueMin = tmpOutputValue;

    itr++;
  }

  outputFile.close();


  string gnuplotCMD;
#ifdef _WIN32
  gnuplotCMD = "set terminal wxt size 1600, 300\n";
#elif __APPLE__
  gnuplotCMD = "set terminal aqua size 1600, 300\n";
#else
  gnuplotCMD = "set terminal qt size 1600, 300 persist\n";
#endif
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
  gnuplotCMD = "set xrange[" + to_string(0) + ":" + to_string(itr - 1) + "]\n";
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
  gnuplotCMD = "set yrange[" + to_string(valueMin - 0.2 * (valueMax
        - valueMin)) + ":" + to_string(valueMax + 0.2 * (valueMax - valueMin))
    + "]\n";
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
  gnuplotCMD = "set grid\n";
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());


  tmpOutputString = " ";
  outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());

  for (int motifSetItr = 0; motifSetItr < (int)motifSetPositions_in.size();
      motifSetItr++) {

    if (motifSetPositions_in[motifSetItr].size() == 2)
      tmpOutputString = ", \"Top Pair Motif Locations\"";
    else
      tmpOutputString = ", \"Top Latent Motif Locations\"";
    outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());
  }

  tmpOutputString = "\n";
  outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());


  tmpOutputString = "\"range/similarity\"";
  outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());

  for (int motifSetItr = 0; motifSetItr < (int)d_in.size(); motifSetItr++) {

          ostringstream streamObj;
          streamObj << fixed << setprecision(6) << ceil(1000000.0
              * d_in[motifSetItr]) / 1000000.0;

    tmpOutputString = ", " + streamObj.str();
    outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());
  }

  tmpOutputString = "\n";
  outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());


  tmpOutputString = "\"window size\"";
  outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());

  for (int motifSetItr = 0; motifSetItr < (int)windowSizes_in.size();
      motifSetItr++) {

    tmpOutputString = ", " + to_string(windowSizes_in[motifSetItr]);
    outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());
  }

  tmpOutputString = "\n";
  outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());


  for (int positionItr = 0; positionItr < maxMSSize; positionItr++) {

    tmpOutputString = "\"position " + to_string(positionItr) + "\"";
    outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());

    for (int motifSetItr = 0; motifSetItr < (int)windowSizes_in.size();
        motifSetItr++) {

      if (positionItr < (int)motifSetPositions_in[motifSetItr].size()) {

        tmpOutputString = ", "
          + to_string(motifSetPositions_in[motifSetItr][positionItr]);
        outputFileMotifSets.write(tmpOutputString.c_str(),
            tmpOutputString.size());

        if (motifSetPositions_in[motifSetItr].size() != 2) {

          gnuplotCMD = "set obj rect from "
            + to_string(motifSetPositions_in[motifSetItr][positionItr])
            + ", graph 0 to "
            + to_string(motifSetPositions_in[motifSetItr][positionItr]
                + windowSizes_in[motifSetItr] - 1)
            + ", graph 1 back fc rgb \"#c5c5c5\" "
            + "fs border rgb \"#c5c5c5\"\n";
          outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
        }
      }
      else {

        tmpOutputString = ", ";
        outputFileMotifSets.write(tmpOutputString.c_str(),
            tmpOutputString.size());
      }
    }

    tmpOutputString = "\n";
    outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());
  }

  for (int positionItr = 0; positionItr < maxMSSize; positionItr++) {

    for (int motifSetItr = 0; motifSetItr < (int)windowSizes_in.size();
        motifSetItr++) {

      if (positionItr < (int)motifSetPositions_in[motifSetItr].size()) {

        if (motifSetPositions_in[motifSetItr].size() == 2) {

          gnuplotCMD = "set obj rect from "
            + to_string(motifSetPositions_in[motifSetItr][positionItr])
            + ", graph 0 to "
            + to_string(motifSetPositions_in[motifSetItr][positionItr]
                + windowSizes_in[motifSetItr] - 1)
            + ", graph 1 back fc rgb \"#656565\" "
            + "fs pattern 1 border rgb \"#656565\"\n";
          outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
        }
      }
    }
  }


  gnuplotCMD = "plot \"" + basicOutputFileName + "_"
    + to_string(outputFolderNumber) + ".csv\" title \"" + timeSeriesName
    + " " + to_string(outputFolderNumber) + "\" with lines";
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());

#ifdef _WIN32
#elif __APPLE__
#else
  gnuplotCMD = "\n\npause mouse close";
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
#endif
  outputFileGNUPlotScript.close();

  outputFileMotifSets.close();
}

void OutputGenerator::setTimeSeriesName(const string &timeSeriesName_in) {

  timeSeriesName = timeSeriesName_in;
}
