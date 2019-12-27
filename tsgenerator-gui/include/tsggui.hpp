///\file tsggui.hpp
///
///\brief Header file of the TSGenerator GUI.
///
///This is the header file of the TSGenerator GUI.

#include <tsgenerator.hpp>
#include <tsgtypes.hpp>
#include <outputgenerator.hpp>
#include <global.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <QObject>
#include <QApplication>
#include <QMainWindow>
#include <QWidget>
#include <QGridLayout>
#include <QtCharts>
#include <QListWidget>
#include <QStringList>
#include <QPushButton>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QRectF>


class TsgGui : public QApplication {

  Q_OBJECT

private:

  //default parameters
  bool running = false;
  int length = 1000;
  int window = 30;
  double delta = 1.0;
  double noise = 2.0;
  tsg::word type = "box";
  int motifSize = 3;
  double height = 10.0;
  double step = 1.0;
  int times = 3;
  double maxi = 20.0;
  tsg::word method = "boundedNormalRandomWalk";
  tsg::word gen = "latent motif";
  int selMotif = 0;

  //final product
  tsg::rseq timeSeries;
  tsg::rseq dVector;
  tsg::iseq windows;
  tsg::iseqs motifPositions;

  //gui components
  QMainWindow gui;
  QWidget pane;
  QGridLayout layout;
  QtCharts::QChartView tsChartView;
  QtCharts::QChart tsChart;
  QtCharts::QLineSeries tsSeries;
  QtCharts::QChartView motifChartView;
  QtCharts::QChart motifChart;
  QtCharts::QLineSeries motifSeries;
  QListWidget motifList;
  QStringList motifs;
  QLabel idxLabel;
  QPushButton startButton;
  QPushButton saveButton;
  QLabel lengthLabel;
  QLineEdit lengthText;
  QLabel windowLabel;
  QLineEdit windowText;
  QLabel sizeLabel;
  QLineEdit sizeText;
  QLabel deltaLabel;
  QLineEdit deltaText;
  QLabel noiseLabel;
  QLineEdit noiseText;
  QLabel heightLabel;
  QLineEdit heightText;
  QLabel stepLabel;
  QLineEdit stepText;
  QLabel timesLabel;
  QLineEdit timesText;
  QLabel maxiLabel;
  QLineEdit maxiText;
  QLabel typeLabel;
  QComboBox typeDrop;
  QStringList types;
  QLabel methodLabel;
  QComboBox methodDrop;
  QStringList methods;
  QLabel genLabel;
  QComboBox genDrop;
  QStringList gens;

public:

  ///\brief The TsgGui constructor.
  ///
  ///\param [in] argc Hands over the number of arguments.
  ///\param [in] *argv[] Hands over the arguments.
  ///
  ///This is the typical main function of C++. The main function is called at
  ///the program start. First all non-local objects with static storage
  ///duration are initialized. We use the main function to parse the argument
  ///list to set up and to start the TSGenerator. Also a configuration file
  ///is created, if a configuration file do not exist.
  TsgGui(int argc, char *argv[]);

  ///\brief The TsgGui destructor.
  ///
  ///Does actually nothing.
  ~TsgGui();

public slots:
  ///\brief This function generates the time series.
  ///
  ///This function generates the time series using the tsgenerator library.
  void generateTS();

  ///\brief This function stores the times series.
  ///
  ///This function stores the ts series and metadata into a seperate folder
  ///using the output generator.
  void saveTS();

  ///\brief This function replots the motif.
  ///
  ///This function replots the selected motif subsequence in the time series.
  void plotMotif();
};

