///\file tsggui.hpp
///
///\brief Header file of the TSGenerator GUI.
///
///This is the header file of the TSGenerator GUI.

#include <tsgenerator.hpp>
#include <tsgtypes.hpp>
#include <outputgenerator.hpp>
#include <global.hpp>
#include <chart.hpp>
#include <chartview.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <QObject>
#include <QApplication>
#include <QMainWindow>
#include <QWidget>
#include <QGridLayout>
#include <QListWidget>
#include <QStringList>
#include <QPushButton>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QRectF>
#include <QMenuBar>
#include <QMenu>
#include <QAction>
#include <QMessageBox>
#include <QLineSeries>
#include <QAreaSeries>
#include <QPen>
#include <QBrush>


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
  QMenuBar menuBar;
  QMenu helpMenu;
  QAction genAct;
  QMessageBox genMes;
  QAction typeAct;
  QMessageBox typeMes;
  QAction methodAct;
  QMessageBox methodMes;
  QAction lengthAct;
  QMessageBox lengthMes;
  QAction windowAct;
  QMessageBox windowMes;
  QAction sizeAct;
  QMessageBox sizeMes;
  QAction noiseAct;
  QMessageBox noiseMes;
  QAction deltaAct;
  QMessageBox deltaMes;
  QAction heightAct;
  QMessageBox heightMes;
  QAction stepAct;
  QMessageBox stepMes;
  QAction timesAct;
  QMessageBox timesMes;
  QAction maxiAct;
  QMessageBox maxiMes;
  QWidget pane;
  QGridLayout layout;
  ChartView tsChartView;
  Chart tsChart;
  QPen tsPen;
  QPen motifPen;
  QBrush brush;
  QLineSeries upperSeries;
  QLineSeries lowerSeries;
  QAreaSeries marker;
  QLineSeries tsSeries;
  ChartView motifChartView;
  Chart motifChart;
  QLineSeries motifSeries;
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

  ///\brief This function opens the generator help.
  ///
  ///This function opens a message box with a help text for the generator
  ///option.
  void showGenHelp();

  ///\brief This function opens the type help.
  ///
  ///This function opens a message box with a help text for the type option.
  void showTypeHelp();

  ///\brief This function opens the method help.
  ///
  ///This function opens a message box with a help text for the method option.
  void showMethodHelp();

  ///\brief This function opens the length help.
  ///
  ///This function opens a message box with a help text for the length option.
  void showLengthHelp();

  ///\brief This function opens the window help.
  ///
  ///This function opens a message box with a help text for the window option.
  void showWindowHelp();

  ///\brief This function opens the size help.
  ///
  ///This function opens a message box with a help text for the size option.
  void showSizeHelp();

  ///\brief This function opens the noise help.
  ///
  ///This function opens a message box with a help text for the noise option.
  void showNoiseHelp();

  ///\brief This function opens the delta help.
  ///
  ///This function opens a message box with a help text for the delta option.
  void showDeltaHelp();

  ///\brief This function opens the height help.
  ///
  ///This function opens a message box with a help text for the height option.
  void showHeightHelp();

  ///\brief This function opens the step help.
  ///
  ///This function opens a message box with a help text for the step option.
  void showStepHelp();

  ///\brief This function opens the times help.
  ///
  ///This function opens a message box with a help text for the times option.
  void showTimesHelp();

  ///\brief This function opens the maxi help.
  ///
  ///This function opens a message box with a help text for the maxi option.
  void showMaxiHelp();
};

