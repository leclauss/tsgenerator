﻿///\file tsggui.hpp
///
///\brief Header file of the TSGenerator GUI.
///
///This is the header file of the TSGenerator GUI.

#include <tsgenerator.hpp>
#include <fstream>
#include <tsgtypes.hpp>
#include <outputgenerator.hpp>
#include <global.hpp>
#include <chart.hpp>
#include <chartview.hpp>
#include <canvas.hpp>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <QObject>
#include <QApplication>
#include <QMainWindow>
#include <QWidget>
#include <QGridLayout>
#include <QList>
#include <QListWidget>
#include <QString>
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
#include <QFileDialog>
#include <QLineSeries>
#include <QAreaSeries>
#include <QPen>
#include <QBrush>
#include <QStringList>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QVBoxLayout>
#include <QDialog>
#include <QDialogButtonBox>


class TsgGui : public QApplication {

  Q_OBJECT

private:

  //default parameters
  int length = tsg::defaultLength;
  int window = tsg::defaultWindow;
  double delta = tsg::defaultDelta;
  double noise = tsg::defaultNoise;
  tsg::word type = tsg::defaultType;
  int motifSize = tsg::defaultMotifSize;
  double height = tsg::defaultHeight;
  double step = tsg::defaultStep;
  int times = tsg::defaultTimes;
  double maxi = tsg::defaultMaxi;
  int smaller = tsg::defaultSmaller;
  tsg::word method = tsg::defaultMethod;
  tsg::word gen = tsg::defaultGen;

  //system variables
  bool running = false;
  std::map<tsg::word, tsg::iseqs> disc;
  tsg::rseq customShape;

  //final product
  tsg::rseq timeSeries;
  tsg::rseq dVector;
  tsg::iseq windows;
  tsg::iseqs motifPositions;
  tsg::rseqs motif;

  //gui components
  QMainWindow gui;
  QMenuBar menuBar;
  QMenu fileMenu;
  QAction openAct;
  QFileDialog openFileDialog;
  QAction saveAct;
  QMenu infoMenu;
  QAction copywriteAct;
  QMessageBox copywriteMes;
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
  QAction smallerAct;
  QMessageBox smallerMes;
  QWidget pane;
  QGridLayout layout;
  QGroupBox tsBox;
  QVBoxLayout tsLayout;
  ChartView tsChartView;
  Chart tsChart;
  QPen tsPen;
  QPen motifPen;
  QBrush brush;
  std::vector<QLineSeries> upperSeriess;
  std::vector<QLineSeries> lowerSeriess;
  std::vector<QAreaSeries> markers;
  QLineSeries tsSeries;
  QGroupBox motifBox;
  QVBoxLayout motifLayout;
  ChartView motifChartView;
  Chart motifChart;
  QGroupBox motifListBox;
  QVBoxLayout motifListLayout;
  std::vector<QLineSeries> motifSeries;
  QListWidget motifList;
  QStringList motifLocs;
  QPushButton startButton;
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
  QLabel smallerLabel;
  QLineEdit smallerText;
  QLabel typeLabel;
  QComboBox typeDrop;
  QStringList types;
  QLabel methodLabel;
  QComboBox methodDrop;
  QStringList methods;
  QLabel genLabel;
  QComboBox genDrop;
  QStringList gens;
  QLabel distLabel;
  QLineEdit distText;
  QLabel rangeLabel;
  QLineEdit rangeText;
  QGroupBox browseBox;
  QHBoxLayout browseLayout;
  QLabel browseLabel;
  QLineEdit browseText;
  QPushButton browseButton;

  //draw dialog box
  QDialog customDia;
  QVBoxLayout customLayout;
  Canvas customCanvas;
  QDialogButtonBox customDiaButtonBox;

public:

  ///\brief The TsgGui constructor.
  ///
  ///\param [in] &argc Hands over the number of arguments.
  ///\param [in] *argv[] Hands over the arguments.
  ///
  ///This is the typical main function of C++. The main function is called at
  ///the program start. First all non-local objects with static storage
  ///duration are initialized. We use the main function to parse the argument
  ///list to set up and to start the TSGenerator. Also a configuration file
  ///is created, if a configuration file do not exist.
  TsgGui(int &argc, char *argv[]);

  ///\brief The TsgGui destructor.
  ///
  ///Does actually nothing.
  ~TsgGui();

  ///\brief This function loads data into the gui.
  ///
  ///This function loads the data into the gui including plotting.
  void loadData();

public slots:
  ///\brief This function generates the time series.
  ///
  ///This function generates the time series using the tsgenerator library.
  void generateTS();

  ///\brief This function opens data.
  ///
  ///This function opens the time series and metadata from a folder
  ///using the output generator.
  void openData();

  ///\brief This function stores data.
  ///
  ///This function stores the time series and metadata into a seperate folder
  ///using the output generator.
  void saveData();

  ///\brief This function replots the motif.
  ///
  ///This function replots the selected motif subsequence in the time series.
  void plotMotif();

  ///\brief This function opens the copywrite.
  ///
  ///This function opens a message box with the copywrite text for the
  ///generator option.
  void showCopywrite();

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

  ///\brief This function opens the smaller option help text.
  ///
  ///This function opens a message box with a help text for the smaller option.
  void showSmallerHelp();

  ///\brief This function opens a drawing window.
  ///
  ///\param [in] index_in Hands over the index of the type drop down list.
  ///
  ///This function opens a box where the user may draw a custom motif shape.
  void typeSelected(const int index_in);

  ///\brief This function creates a custom shape from a drawing.
  ///
  ///This function creates a custom shape from the drawing of the user.
  void createCustomShape();

  ///\brief This function closes the drawing window.
  ///
  ///This function closes the box where the user may draw a custom motif shape.
  void closeCustomDia();
};
