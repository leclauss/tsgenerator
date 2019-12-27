///\file outputgenerator.cpp
///
///\brief File contains the output generator class definition.
///
///This is the source file of the output generator. The output generator writes
///the time series and the motif set positions into output files.

#include <outputgenerator.hpp>

OutputGenerator::OutputGenerator() : OutputGenerator("time_series") {}

OutputGenerator::OutputGenerator(const std::string &outputFileName_in)
  : outputFileName("/" + outputFileName_in) {

  std::string outputDirName = "./time_series_";

  while (std::filesystem::exists(outputDirName
        + std::to_string(outputFolderNumber)))
    outputFolderNumber++;

  std::filesystem::create_directory(outputDirName
      + std::to_string(outputFolderNumber));

  basicOutputFileName = outputFileName;
  outputFileMotifSetsName = outputDirName + std::to_string(outputFolderNumber)
    + outputFileName + "_meta";
  outputFileGNUPlotScriptName = outputDirName
    + std::to_string(outputFolderNumber) + outputFileName + "_plot";
  outputFileName = outputDirName + std::to_string(outputFolderNumber)
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

  std::string tmpTimeSeriesFileName = outputFileName + "_"
    + std::to_string(outputFolderNumber) + ".csv";
  std::string tmpFileMotifSetPositionsName = outputFileMotifSetsName + "_"
    + std::to_string(outputFolderNumber) + ".csv";
  std::string tmpGNUPlotScriptName = outputFileGNUPlotScriptName + "_"
    + std::to_string(outputFolderNumber) + ".plt";

  if (std::filesystem::exists(tmpTimeSeriesFileName)) {

    std::cerr << "ERROR: Time series file already exist!" << std::endl;
    throw(EXIT_FAILURE);
  }

  if (std::filesystem::exists(tmpFileMotifSetPositionsName)) {

    std::cerr << "ERROR: Motif positions file already exist!" << std::endl;
    throw(EXIT_FAILURE);
  }

  if (std::filesystem::exists(tmpGNUPlotScriptName)) {

    std::cerr << "ERROR: GNUPlot script file already exist!" << std::endl;
    throw(EXIT_FAILURE);
  }

  outputFile.open(tmpTimeSeriesFileName, std::ios::out);
  outputFileMotifSets.open(tmpFileMotifSetPositionsName, std::ios::out);
  outputFileGNUPlotScript.open(tmpGNUPlotScriptName, std::ios::out);
}

void OutputGenerator::printTimeSeriesHorizontal(const std::vector<double>
    &timeSeries_in, const std::vector<double> &d_in, const std::vector<int>
    &windowSizes_in, const std::vector<std::vector<int>> &motifSetPositions_in)
{

  if (timeSeries_in.size() == 0) {

    std::cerr << "ERROR: Time series is empty." << std::endl;
    throw(EXIT_FAILURE);
  }

  if (d_in.size() == 0) {

    std::cerr << "ERROR: Distance vector is empty." << std::endl;
    throw(EXIT_FAILURE);
  }

  if (d_in.size() != motifSetPositions_in.size()) {

    std::cerr << "ERROR: Distance vector and motif set positions vector" <<
      " have different size." << std::endl;
    throw(EXIT_FAILURE);
  }

  if (windowSizes_in.size() == 0) {

    std::cerr << "ERROR: Window size vector is empty." << std::endl;
    throw(EXIT_FAILURE);
  }

  if (windowSizes_in.size() != d_in.size()) {

    std::cerr << "ERROR: Window size vector and distance vector have" <<
      " different size." << std::endl;
    throw(EXIT_FAILURE);
  }

  //check if motif set positions vector is empty
  int maxMSSize = 0;

  for (int motifSetItr = 0; motifSetItr < (int)d_in.size(); motifSetItr++) {

    if ((motifSetPositions_in[motifSetItr]).size() == 0) {

      std::cerr << "ERROR: Motif set positions vector is empty." << std::endl;
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
  std::string tmpOutputString;

  tmpOutputValue = timeSeries_in[0];
  tmpOutputString = std::to_string(tmpOutputValue);
  outputFile.write(tmpOutputString.c_str(), tmpOutputString.size());

  while (itr < (int)timeSeries_in.size()) {

    tmpOutputValue = timeSeries_in[itr];
    tmpOutputString = ", " + std::to_string(tmpOutputValue);
    outputFile.write(tmpOutputString.c_str(), tmpOutputString.size());

    if (tmpOutputValue > valueMax)
      valueMax = tmpOutputValue;
    if (tmpOutputValue < valueMin)
      valueMin = tmpOutputValue;

    itr++;
  }

  outputFile.close();


  tmpOutputString = "\n";

  std::string gnuplotCMD;
#ifdef _WIN32
  gnuplotCMD = "set terminal wxt size 1600, 300\n";
#elif __APPLE__
  gnuplotCMD = "set terminal aqua size 1600, 300\n";
#else
  gnuplotCMD = "set terminal qt size 1600, 300 persist\n";
#endif
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
  gnuplotCMD = "set xrange[" + std::to_string(0) + ":" + std::to_string(itr
      - 1) + "]\n";
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
  gnuplotCMD = "set yrange[" + std::to_string(valueMin - 0.2 * (valueMax
        - valueMin)) + ":" + std::to_string(valueMax + 0.2 * (valueMax
          - valueMin)) + "]\n";
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
  gnuplotCMD = "set grid\n";
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());


  int maxMotifSetSize = 0;

  for (std::vector<int> item : motifSetPositions_in)
    if ((int)item.size() > maxMotifSetSize)
      maxMotifSetSize = item.size();

  tmpOutputString = " , \"range/similarity\"";
  outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());

  tmpOutputString = ", \"window size\"";
  outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());

  for (int vectorItr = 0; vectorItr < maxMotifSetSize; vectorItr++) {

    tmpOutputString = ", \"position " + std::to_string(vectorItr) + "\"";
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

          std::ostringstream streamObj;
          streamObj << std::fixed << std::setprecision(6) <<
            std::ceil(1000000.0 * d_in[vectorItr]) / 1000000.0;

    tmpOutputString = ", " + streamObj.str();
    outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());

    tmpOutputString =  ", " + std::to_string(windowSizes_in[vectorItr]);
    outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());

    for (int positionItr = 0; positionItr
        < (int)(motifSetPositions_in[vectorItr]).size(); positionItr++) {

      tmpOutputString = ", "
        + std::to_string((motifSetPositions_in[vectorItr])[positionItr]);
      outputFileMotifSets.write(tmpOutputString.c_str(),
          tmpOutputString.size());

      if (motifSetPositions_in[vectorItr].size() != 2) {

        gnuplotCMD = "set obj rect from "
          + std::to_string(motifSetPositions_in[vectorItr][positionItr])
          + ", graph 0 to "
          + std::to_string(motifSetPositions_in[vectorItr][positionItr]
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
          + std::to_string(motifSetPositions_in[vectorItr][positionItr])
          + ", graph 0 to "
          + std::to_string(motifSetPositions_in[vectorItr][positionItr]
              + windowSizes_in[vectorItr] - 1)
          + ", graph 1 back fc rgb \"#656565\" "
          + "fs pattern 1 border rgb \"#656565\"\n";
        outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
      }
    }
  }


  gnuplotCMD = "plot \"" + basicOutputFileName + "_"
    + std::to_string(outputFolderNumber) + ".csv\" matrix title \""
    + timeSeriesName + " " + std::to_string(outputFolderNumber) + "\" with "
    "lines";
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

void OutputGenerator::printTimeSeriesVertical(const std::vector<double>
    &timeSeries_in, const std::vector<double> &d_in, const std::vector<int>
    &windowSizes_in, const std::vector<std::vector<int>> &motifSetPositions_in) {

  if (timeSeries_in.size() == 0) {

    std::cerr << "ERROR: Time series is empty." << std::endl;
    throw(EXIT_FAILURE);
  }

  if (d_in.size() == 0) {

    std::cerr << "ERROR: Distance vector is empty." << std::endl;
    throw(EXIT_FAILURE);
  }

  if (d_in.size() != motifSetPositions_in.size()) {

    std::cerr << "ERROR: Distance vector and motif set positions vector" <<
      " have different size." << std::endl;
    throw(EXIT_FAILURE);
  }

  if (windowSizes_in.size() == 0) {

    std::cerr << "ERROR: Window size vector is empty." << std::endl;
    throw(EXIT_FAILURE);
  }

  if (windowSizes_in.size() != d_in.size()) {

    std::cerr << "ERROR: Window size vector and distance vector have" <<
      " different size." << std::endl;
    throw(EXIT_FAILURE);
  }

  //check if there is an empty motif set vector
  int maxMSSize = 0;

  for (int motifSetItr = 0; motifSetItr < (int)d_in.size(); motifSetItr++) {

    if ((int)(motifSetPositions_in[motifSetItr]).size() == 0) {

      std::cerr << "ERROR: Motif set positions vector is empty." << std::endl;
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
  std::string tmpOutputString;

  while (itr < (int)timeSeries_in.size()) {

    tmpOutputValue = timeSeries_in[itr];
    tmpOutputString = std::to_string(tmpOutputValue) + "\n";
    outputFile.write(tmpOutputString.c_str(), tmpOutputString.size());

    if (tmpOutputValue > valueMax)
      valueMax = tmpOutputValue;
    if (tmpOutputValue < valueMin)
      valueMin = tmpOutputValue;

    itr++;
  }

  outputFile.close();


  std::string gnuplotCMD;
#ifdef _WIN32
  gnuplotCMD = "set terminal wxt size 1600, 300\n";
#elif __APPLE__
  gnuplotCMD = "set terminal aqua size 1600, 300\n";
#else
  gnuplotCMD = "set terminal qt size 1600, 300 persist\n";
#endif
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
  gnuplotCMD = "set xrange[" + std::to_string(0) + ":" + std::to_string(itr
      - 1) + "]\n";
  outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
  gnuplotCMD = "set yrange[" + std::to_string(valueMin - 0.2 * (valueMax
        - valueMin)) + ":" + std::to_string(valueMax + 0.2 * (valueMax
          - valueMin)) + "]\n";
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

          std::ostringstream streamObj;
          streamObj << std::fixed << std::setprecision(6) <<
            std::ceil(1000000.0 * d_in[motifSetItr]) / 1000000.0;

    tmpOutputString = ", " + streamObj.str();
    outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());
  }

  tmpOutputString = "\n";
  outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());


  tmpOutputString = "\"window size\"";
  outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());

  for (int motifSetItr = 0; motifSetItr < (int)windowSizes_in.size();
      motifSetItr++) {

    tmpOutputString = ", " + std::to_string(windowSizes_in[motifSetItr]);
    outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());
  }

  tmpOutputString = "\n";
  outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());


  for (int positionItr = 0; positionItr < maxMSSize; positionItr++) {

    tmpOutputString = "\"position " + std::to_string(positionItr) + "\"";
    outputFileMotifSets.write(tmpOutputString.c_str(), tmpOutputString.size());

    for (int motifSetItr = 0; motifSetItr < (int)windowSizes_in.size();
        motifSetItr++) {

      if (positionItr < (int)motifSetPositions_in[motifSetItr].size()) {

        tmpOutputString = ", "
          + std::to_string(motifSetPositions_in[motifSetItr][positionItr]);
        outputFileMotifSets.write(tmpOutputString.c_str(),
            tmpOutputString.size());

        if (motifSetPositions_in[motifSetItr].size() != 2) {

          gnuplotCMD = "set obj rect from "
            + std::to_string(motifSetPositions_in[motifSetItr][positionItr])
            + ", graph 0 to "
            + std::to_string(motifSetPositions_in[motifSetItr][positionItr]
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
            + std::to_string(motifSetPositions_in[motifSetItr][positionItr])
            + ", graph 0 to "
            + std::to_string(motifSetPositions_in[motifSetItr][positionItr]
                + windowSizes_in[motifSetItr] - 1)
            + ", graph 1 back fc rgb \"#656565\" "
            + "fs pattern 1 border rgb \"#656565\"\n";
          outputFileGNUPlotScript.write(gnuplotCMD.c_str(), gnuplotCMD.size());
        }
      }
    }
  }


  gnuplotCMD = "plot \"" + basicOutputFileName + "_"
    + std::to_string(outputFolderNumber) + ".csv\" title \"" + timeSeriesName
    + " " + std::to_string(outputFolderNumber) + "\" with lines";
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

void OutputGenerator::setTimeSeriesName(const std::string &timeSeriesName_in) {

  timeSeriesName = timeSeriesName_in;
}
