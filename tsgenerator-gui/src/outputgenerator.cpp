///\file outputgenerator.cpp
///
///\brief File contains the output generator class definition.
///
///This is the source file of the output generator. The output generator writes
///the time series and the motif set positions into output files.

#include <outputgenerator.hpp>


OutputGenerator::OutputGenerator(const std::string fileName_in, const
    std::string path_in, const std::string delimiter_in)
  : delimiter(delimiter_in) {

  dir = "time_series_";

  if (path_in.back() == '/')
    dir = path_in + dir;
  else
    dir = path_in + "/" + dir;

  baseFileName = fileName_in;
}

OutputGenerator::~OutputGenerator() {

  close();
}

void OutputGenerator::open() {

  close();

  while (std::filesystem::exists(dir + std::to_string(folderNumber)))
    folderNumber++;

  std::filesystem::create_directory(dir + std::to_string(folderNumber));

  metaFileName = dir + std::to_string(folderNumber) + "/" + baseFileName
    + "_meta";
  tsFileName = dir + std::to_string(folderNumber) + "/" + baseFileName;

  tsFileName += "_" + std::to_string(folderNumber) + ".csv";
  metaFileName += "_" + std::to_string(folderNumber) + ".csv";

  if (std::filesystem::exists(tsFileName)) {

    std::cerr << "ERROR: Time series file already exist!" << std::endl;
    throw(EXIT_FAILURE);
  }

  if (std::filesystem::exists(metaFileName)) {

    std::cerr << "ERROR: Motif positions file already exist!" << std::endl;
    throw(EXIT_FAILURE);
  }

  tsFile.open(tsFileName, std::ios::out);
  metaFile.open(metaFileName, std::ios::out);
}

void OutputGenerator::close() {

  if (tsFile.is_open())
    tsFile.close();

  if (metaFile.is_open())
    metaFile.close();
}

void OutputGenerator::printTimeSeriesHorizontal(const std::vector<double>
    &timeSeries_in) {

  if (!timeSeries_in.empty())
    tsFile << timeSeries_in[0];

  for (int i = 1; i < (int)timeSeries_in.size(); i++)
    tsFile << delimiter << timeSeries_in[i];
}

void OutputGenerator::printTimeSeriesVertical(const std::vector<double>
    &timeSeries_in) {

  for (int i = 0; i < (int)timeSeries_in.size(); i++)
    tsFile << timeSeries_in[i] << std::endl;
}

void OutputGenerator::printMetaLine(const std::vector<int> &data_in, const
    std::string name_in) {

  metaFile << name_in;

  for (int i = 0; i < (int)data_in.size(); i++)
    metaFile << delimiter << data_in[i];

  metaFile << std::endl;
}

void OutputGenerator::printMetaLine(const std::vector<double> &data_in, const
    std::string name_in) {

  metaFile << name_in;

  for (int i = 0; i < (int)data_in.size(); i++)
    metaFile << delimiter << std::fixed << std::setprecision(6)
      << std::ceil(1000000.0 * data_in[i]) / 1000000.0;

  metaFile << std::endl;
}

void OutputGenerator::printMetaLine(const std::vector<std::string> &data_in,
    const std::string name_in) {

  metaFile << name_in;

  for (int i = 0; i < (int)data_in.size(); i++)
    metaFile << delimiter << data_in[i];

  metaFile << std::endl;
}

void OutputGenerator::setFileName(const std::string &fileName_in) {

  baseFileName = "/" + fileName_in;
  metaFileName = "./time_series_" + std::to_string(folderNumber)
    + "/" + fileName_in + "_meta";
  tsFileName = "./time_series_" + std::to_string(folderNumber)
    + "/" + fileName_in;
}
