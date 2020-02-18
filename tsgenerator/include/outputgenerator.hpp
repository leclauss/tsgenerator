///\file outputgenerator.hpp
///
///\brief File contains the output generator class declaration.
///
///This is the header file of the output generator. The output generator writes
///the time series and the motif positions into output files.

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
#include <filesystem>


///\brief This class represents a output generator.
///
///The OutputGenerator class offers a horizontal and vertical data print
///functions. Therefore, The time series and the motif positions list will be
///printed either horizontal or vertical into the output files.
class OutputGenerator {

protected:

  ///\brief Stores the time series output file pointer.
  ///
  ///This variable stores the pointer pointing at the time series output file.
  std::ofstream tsFile;

  ///\brief Stores the meta output file pointer.
  ///
  ///This variable stores the pointer pointing at the output file containing
  ///the meta data of the time series.
  std::ofstream metaFile;

  ///\brief Stores the path to the output folder.
  ///
  ///This variable stores the path to the output folder used to store the
  ///output files.
  std::string dir;

  ///\brief Stores the name of the meta output file.
  ///
  ///This variable stores the name of the file containing the meta data of the
  ///time series.
  std::string metaFileName;

  ///\brief Stores the prefix of the time series file.
  ///
  ///This variable stores the prefix of the time series output file.
  std::string tsFileName;

  ///\brief Stores the time series file name.
  ///
  ///This variable stores the time series output file name.
  std::string baseFileName;

  ///\brief Stores the number of the output file.
  ///
  ///This variable stores the nubmer of the current output file.
  int folderNumber = 0;

  ///\brief Stores the delimiter sequence.
  ///
  ///This variable stores the delimiter symbol dividing the values in the
  ///output files.
  std::string delimiter;

public:

  ///\brief Initializes the output generator.
  ///
  ///\param [in] delimiter_in Sets the delimiter sequence.
  ///
  ///The constructor creates a new output directory path and sets up the file
  ///name with default file name "time_series". It chooses another directory
  ///path, if a directory with the same name already exists.
  OutputGenerator(const std::string delimiter_in = ",");

  ///\brief Initializes the output generator, i.e. creates output directory.
  ///
  ///\param [in] outputFileName_in Hands over the output file name.
  ///\param [in] delimiter_in Hands over the delimiter sequence.
  ///\param [in] path_in Hands over the path for the output files.
  ///
  ///The constructor creates a new output directory path and sets up the file
  ///name with default file name "time_series". It chooses another directory
  ///path, if a directory with the same name already exists.
  OutputGenerator(const std::string outputFileName_in, const std::string
      delimiter_in = ",", const std::string path_in = ".");

  ///\brief Destroys the file pointers.
  ///
  ///The destructor destroys the time series file pointer and the motif set
  ///positions file pointer.
  ~OutputGenerator();

  ///\brief Opens files for the output of the time series and meta data.
  ///
  ///This function opens the files for the time series and the meta data
  ///output.
  void open();

  ///\brief Closes files for the output of the time series and meta data.
  ///
  ///This function closes the files for the time series and meta data output.
  void close();

  ///\brief Writes the data horizontal into the output file.
  ///
  ///\param [in] &timeSeries_in Hands over the time series that output
  ///generator writes into output file.
  ///
  ///This function writes the time series into the output file divided by the
  ///delimiter.
  void printTimeSeriesHorizontal(const std::vector<double> &timeSeries_in);

  ///\brief Writes the data vertical into the output files.
  ///
  ///\param [in] &timeSeries_in Hands over the time series that output
  ///generator writes into output file.
  ///
  ///This function writes the time series values into the output file line by
  ///line.
  void printTimeSeriesVertical(const std::vector<double> &timeSeries_in);

  ///\brief Writes a line into the meta data output file.
  ///
  ///\param [in] &data_in Hands over the values to output.
  ///\param [in] name_in Hands over the name of the data line.
  ///
  ///The print meta line function writes a line name_in delimiter data_in[0]
  ///data_in[1] ... into the meta data output file.
  void printMetaLine(const std::vector<int> &data_in, const std::string
      name_in);

  ///\brief Writes a line into the meta data output file.
  ///
  ///\param [in] &data_in Hands over the values to output.
  ///\param [in] name_in Hands over the name of the data line.
  ///
  ///The print meta line function writes a line name_in delimiter data_in[0]
  ///data_in[1] ... into the meta data output file. The real values are rounded
  ///to 6 digits after the decimal point.
  void printMetaLine(const std::vector<double> &data_in, const std::string
      name_in);

  ///\brief Writes a line into the meta data output file.
  ///
  ///\param [in] &data_in Hands over the values to output.
  ///\param [in] name_in Hands over the name of the data line.
  ///
  ///The print meta line function writes a line name_in delimiter data_in[0]
  ///data_in[1] ... into the meta data output file.
  void printMetaLine(const std::vector<std::string> &data_in, const std::string
      name_in);

  ///\brief This is the file name setter function.
  ///
  ///\param [in] &fileName_in Hands over the file name.
  ///
  ///This function sets the name of the output files.
  void setFileName(const std::string &fileName_in);
};

#endif
