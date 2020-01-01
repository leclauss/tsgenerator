///\file tsggui.cpp
///
///\brief Source file of the TSGenerator GUI.
///
///This is the source file of the TSGenerator GUI.

#include <tsggui.hpp>


TsgGui::TsgGui(int argc, char *argv[]) : QApplication(argc, argv) {

  //help message box

  //set initialize widgets
  helpMenu.setTitle("Help");
  genAct.setText("generator");
  genMes.setWindowTitle("generator option");
  genMes.setText("The generator option specifies the method for injecting "
      "sequences into the time series matching the synthetic motif.");
  typeAct.setText("type");
  typeMes.setWindowTitle("type option");
  typeMes.setText("The type option specifies the shape of the injected motif "
      "and the inserted sequences.");
  methodAct.setText("method");
  methodMes.setWindowTitle("method option");
  methodMes.setText("The method option sets the method for generating the "
      "base time series.");
  lengthAct.setText("length");
  lengthMes.setWindowTitle("length option");
  lengthMes.setText("The length option specifies the length of the time "
      "series. The length is the number of values in the time series.");
  windowAct.setText("window");
  windowMes.setWindowTitle("window option");
  windowMes.setText("The window option sets the sequence size of the "
      "sequences injected into the time series non-self matching the motif.");
  sizeAct.setText("size");
  sizeMes.setWindowTitle("size option");
  sizeMes.setText("The size option sets the number of subsequences matching "
      "the motif sequence.");
  noiseAct.setText("noise");
  noiseMes.setWindowTitle("noise option");
  noiseMes.setText("The noise option option is the value added to the base "
      "time series and the injected sequences.");
  deltaAct.setText("delta");
  deltaMes.setWindowTitle("delta option");
  deltaMes.setText("The delta option is the difference between two "
      "consecutive values in the time series.");
  heightAct.setText("height");
  heightMes.setWindowTitle("height option");
  heightMes.setText("The height option sets the height of the motif, i.e. the "
      "maximum difference between two values of the base motif.");
  stepAct.setText("step");
  stepMes.setWindowTitle("step option");
  stepMes.setText("The step option sets the maximum step size from two "
      "consecutive values when creating a splined base time series.");
  timesAct.setText("times");
  timesMes.setWindowTitle("times option");
  timesMes.setText("The times option sets the number of values of the "
      "repeating pattern when a splined base time series is generated");
  maxiAct.setText("maxi");
  maxiMes.setWindowTitle("maxi option");
  maxiMes.setText("The maxi option sets the maximum absolute value of the "
      "base time series.");
  idxLabel.setText("Motif Locations");
  startButton.setText("start");
  saveButton.setText("save");
  lengthLabel.setText("length");
  lengthText.setPlaceholderText(std::to_string(length).c_str());
  windowLabel.setText("window");
  windowText.setPlaceholderText(std::to_string(window).c_str());
  sizeLabel.setText("size");
  sizeText.setPlaceholderText(std::to_string(motifSize).c_str());
  deltaLabel.setText("delta");
  deltaText.setPlaceholderText(std::to_string(delta).c_str());
  noiseLabel.setText("noise");
  noiseText.setPlaceholderText(std::to_string(noise).c_str());
  heightLabel.setText("height");
  heightText.setPlaceholderText(std::to_string(height).c_str());
  stepLabel.setText("step");
  stepText.setPlaceholderText(std::to_string(step).c_str());
  timesLabel.setText("times");
  timesText.setPlaceholderText(std::to_string(times).c_str());
  maxiLabel.setText("maxi");
  maxiText.setPlaceholderText(std::to_string(maxi).c_str());
  typeLabel.setText("type");
  methodLabel.setText("method");
  genLabel.setText("generator");

  //connect the message boxes to the actions
  connect(&genAct, SIGNAL(triggered()), this, SLOT(showGenHelp()));
  connect(&typeAct, SIGNAL(triggered()), this, SLOT(showTypeHelp()));
  connect(&methodAct, SIGNAL(triggered()), this, SLOT(showMethodHelp()));
  connect(&lengthAct, SIGNAL(triggered()), this, SLOT(showLengthHelp()));
  connect(&windowAct, SIGNAL(triggered()), this, SLOT(showWindowHelp()));
  connect(&sizeAct, SIGNAL(triggered()), this, SLOT(showSizeHelp()));
  connect(&noiseAct, SIGNAL(triggered()), this, SLOT(showNoiseHelp()));
  connect(&deltaAct, SIGNAL(triggered()), this, SLOT(showDeltaHelp()));
  connect(&heightAct, SIGNAL(triggered()), this, SLOT(showHeightHelp()));
  connect(&stepAct, SIGNAL(triggered()), this, SLOT(showStepHelp()));
  connect(&timesAct, SIGNAL(triggered()), this, SLOT(showTimesHelp()));
  connect(&maxiAct, SIGNAL(triggered()), this, SLOT(showMaxiHelp()));

  //create menus
  helpMenu.addAction(&genAct);
  helpMenu.addAction(&typeAct);
  helpMenu.addAction(&methodAct);
  helpMenu.addAction(&lengthAct);
  helpMenu.addAction(&windowAct);
  helpMenu.addAction(&sizeAct);
  helpMenu.addAction(&noiseAct);
  helpMenu.addAction(&deltaAct);
  helpMenu.addAction(&heightAct);
  helpMenu.addAction(&stepAct);
  helpMenu.addAction(&timesAct);
  helpMenu.addAction(&maxiAct);
  menuBar.addMenu(&helpMenu);

  gui.setMenuBar(&menuBar);

  //add the types
  types.append("box");
  types.append("triangle");
  types.append("semicircle");
  types.append("trapezoid");
  types.append("positiveflank");
  types.append("negativeflank");
  types.append("sine");
  types.append("cosine");

  typeDrop.addItems(types);

  //add the methods
  methods.append("simpleRandomWalk");
  methods.append("realRandomWalk");
  methods.append("normalRandomWalk");
  methods.append("linearRandomWalk");
  methods.append("boundedSimpleRandomWalk");
  methods.append("boundedRealRandomWalk");
  methods.append("boundedNormalRandomWalk");
  methods.append("boundedLinearRandomWalk");
  methods.append("uniformRandom");
  methods.append("normalRandom");
  methods.append("piecewiseLinearRandom");
  methods.append("splineRepeated");

  methodDrop.addItems(methods);
  methodDrop.setCurrentIndex(6);

  //add the generators
  gens.append("pair motif");
  gens.append("set motif");
  gens.append("latent motif");

  genDrop.addItems(gens);
  genDrop.setCurrentIndex(2);

  //setup buttons
  connect(&startButton, SIGNAL(clicked()), this, SLOT(generateTS()));
  connect(&saveButton, SIGNAL(clicked()), this, SLOT(saveTS()));

  //set subsqeuence marker
  motifPen.setColor(0xffbb00);
  brush.setColor(0xffddaa);
  brush.setStyle(Qt::SolidPattern);
  marker.setPen(motifPen);
  marker.setBrush(brush);
  marker.setUpperSeries(&upperSeries);
  marker.setLowerSeries(&lowerSeries);

  //set time series color
  tsPen.setColor(0x0076ff);
  tsSeries.setPen(tsPen);

  //set motif plot color
  motifSeries.setPen(motifPen);

  //add time series plot
  tsChart.addSeries(&marker);
  tsChart.addSeries(&tsSeries);
  tsChart.createDefaultAxes();
  tsChart.legend()->hide();

  tsChartView.setChart(&tsChart);
  tsChartView.setRenderHint(QPainter::Antialiasing);
  tsChartView.setMinimumSize(1200, 200);
  tsChartView.setMaximumSize(30000, 500);

  //add motif location plot
  motifChart.addSeries(&motifSeries);
  motifChart.createDefaultAxes();
  motifChart.legend()->hide();

  motifChartView.setChart(&motifChart);
  motifChartView.setRenderHint(QPainter::Antialiasing);
  motifChartView.setMinimumSize(20, 250);
  motifChartView.setMaximumSize(30000, 600);

  //add plot motif location trigger
  connect(&motifList, SIGNAL(itemSelectionChanged()), this, SLOT(plotMotif()));

  //layout all items
  layout.addWidget(&genLabel, 0, 0);
  layout.addWidget(&genDrop, 1, 0);
  layout.addWidget(&typeLabel, 0, 1);
  layout.addWidget(&typeDrop, 1, 1);
  layout.addWidget(&methodLabel, 0, 2);
  layout.addWidget(&methodDrop, 1, 2);
  layout.addWidget(&lengthLabel, 0, 3);
  layout.addWidget(&lengthText, 1, 3);
  layout.addWidget(&windowLabel, 0, 4);
  layout.addWidget(&windowText, 1, 4);
  layout.addWidget(&sizeLabel, 0, 5);
  layout.addWidget(&sizeText, 1, 5);
  layout.addWidget(&noiseLabel, 2, 0);
  layout.addWidget(&noiseText, 3, 0);
  layout.addWidget(&deltaLabel, 2, 1);
  layout.addWidget(&deltaText, 3, 1);
  layout.addWidget(&heightLabel, 2, 2);
  layout.addWidget(&heightText, 3, 2);
  layout.addWidget(&stepLabel, 2, 3);
  layout.addWidget(&stepText, 3, 3);
  layout.addWidget(&timesLabel, 2, 4);
  layout.addWidget(&timesText, 3, 4);
  layout.addWidget(&maxiLabel, 2, 5);
  layout.addWidget(&maxiText, 3, 5);
  layout.addWidget(&startButton, 4, 4);
  layout.addWidget(&saveButton, 4, 5);
  layout.addWidget(&tsChartView, 5, 0, 1, -1);
  layout.addWidget(&motifChartView, 6, 0, 2, 3);
  layout.addWidget(&idxLabel, 6, 3);
  layout.addWidget(&motifList, 7, 3);

  pane.setLayout(&layout);

  //build window
  gui.setCentralWidget(&pane);
  gui.resize(800, 300);
  gui.setWindowTitle("tsgenerator");
  gui.grabGesture(Qt::PanGesture);
  gui.grabGesture(Qt::PinchGesture);
  gui.show();
}

TsgGui::~TsgGui() { }

void TsgGui::generateTS() {

  if (!running) {

    running = true;
    bool first = false;
    bool success = false;
    double min, max;
    std::string in;

    //read the input
    in = lengthText.text().toStdString();
    if (in.length() > 0 && in.find_first_not_of("0123456789") ==
        std::string::npos)
      length = std::stoi(in);

    in = windowText.text().toStdString();
    if (in.length() > 0 && in.find_first_not_of("0123456789") ==
        std::string::npos)
      window = std::stoi(in);

    in = deltaText.text().toStdString();
    if (check_if_float(in))
      delta = std::stod(in);

    in = noiseText.text().toStdString();
    if (check_if_float(in))
      noise = std::stod(in);

    in = typeDrop.currentText().toStdString();
    if (in.length())
      type = in;

    in = sizeText.text().toStdString();
    if (in.length() > 0 && in.find_first_not_of("0123456789") ==
        std::string::npos)
      motifSize = std::stoi(in);

    in = heightText.text().toStdString();
    if (check_if_float(in))
      height = std::stod(in);

    in = stepText.text().toStdString();
    if (in.length() > 0 && in.find_first_not_of("0123456789") ==
        std::string::npos)
      step = std::stoi(in);

    in = timesText.text().toStdString();
    if (in.length() > 0 && in.find_first_not_of("0123456789") ==
        std::string::npos)
      times = std::stoi(in);

    in = methodDrop.currentText().toStdString();
    if (in.length())
      method = in;

    in = genDrop.currentText().toStdString();
    if (in.length())
      gen = in;

    in = maxiText.text().toStdString();
    if (check_if_float(in))
      maxi = std::stod(in);

    try {

      //generate the time series
      tsg::TSGenerator tSGenerator(length, window, delta, noise, type,
          motifSize, height, step, times, method, maxi, gen);
      tSGenerator.run(timeSeries, motif, dVector, motifPositions);

      success = true;
    }
    catch (int e) {

      if (e != EXIT_FAILURE)
        std::cerr << "ERROR: Something unexpected happened!" << std::endl;
    }
    catch (...) {

      std::cerr << "ERROR: Something unexpected happened!" << std::endl;
    }

    if (success) {

      //check if first time series
      if (motifs.empty())
        first = true;
      else
        selMotif = 0;

      //highlight motif
      if (!first) {

        upperSeries.clear();
        lowerSeries.clear();
      }

      upperSeries.append(motifPositions[0][selMotif],
          100000000.0);
      upperSeries.append(motifPositions[0][selMotif] + window - 1,
          100000000.0);
      lowerSeries.append(motifPositions[0][selMotif],
          -100000000.0);
      lowerSeries.append(motifPositions[0][selMotif] + window - 1,
          -100000000.0);


      //print ts
      if (!first)
        tsSeries.clear();

      tsSeries.append(0, timeSeries[0]);
      min = timeSeries[0];
      max = min;

      for (int i = 1; i < length; i++) {

        if (timeSeries[i] < min)
          min = timeSeries[i];
        if (timeSeries[i] > max)
          max = timeSeries[i];

        tsSeries.append(i, timeSeries[i]);
      }

      tsChart.axes(Qt::Horizontal)[0]->setRange(0.0, length);
      tsChart.axes(Qt::Vertical)[0]->setRange(min - 0.1 * (max - min), max
          + 0.1 * (max - min));

      //print motif location plot
      if (!first)
        motifSeries.clear();

      motifSeries.append(motifPositions[0][selMotif],
          timeSeries[motifPositions[0][selMotif]]);
      min = timeSeries[motifPositions[0][selMotif]];
      max = min;

      for (int i = motifPositions[0][selMotif] + 1;
          i < motifPositions[0][selMotif] + window; i++) {

        if (timeSeries[i] < min)
          min = timeSeries[i];
        if (timeSeries[i] > max)
          max = timeSeries[i];
        motifSeries.append(i, timeSeries[i]);
      }

      motifChart.axes(Qt::Horizontal)[0]->setRange(motifPositions[0][selMotif],
          motifPositions[0][selMotif] + window - 1);
      motifChart.axes(Qt::Vertical)[0]->setRange(min - 0.1 * (max - min), max
          + 0.1 * (max - min));

      //print motif select list
      if (!first) {

        motifs.clear();
        motifList.clear();
      }

      for (int i = 0; i < (int)motifPositions[0].size(); i++) {

        std::stringstream idx;
        idx << "idx: " << motifPositions[0][i];
        motifs.append(idx.str().c_str());
      }

      motifList.addItems(motifs);
      motifList.setCurrentRow(0);

      selMotif = 0;
    }
    else {
      std::cerr << "ERROR: Couldn't create time series." << std::endl;
    }

    running = false;
  }
}

void TsgGui::saveTS() {

  if (!running) {

    running = true;

    try {

      //output stuff
      OutputGenerator outputFile;

      outputFile.open();
      outputFile.printTimeSeriesVertical(timeSeries);

      if (gen == "set motif") {

        outputFile.printMetaLine((std::vector<int>){ window }, "window");
        outputFile.printMetaLine((std::vector<int>){ motifPositions[0][0] },
            "set motif");
        outputFile.printMetaLine(motifPositions[0], "matchings");
        outputFile.printMetaLine((std::vector<double>){ dVector[0] }, "range");
        outputFile.printMetaLine(motifPositions[1], "pair motif");
        outputFile.printMetaLine((std::vector<double>){ dVector[1] }, "pair "
            "motif distance");
      }
      else if (gen == "latent motif") {

        outputFile.printMetaLine((std::vector<int>){ window }, "window");
        //tbd
        outputFile.printMetaLine((std::vector<double>){}, "latent motif");
        outputFile.printMetaLine(motifPositions[0], "matchings");
        outputFile.printMetaLine((std::vector<double>){ dVector[0] }, "range");
        outputFile.printMetaLine(motifPositions[1], "pair motif");
        outputFile.printMetaLine((std::vector<double>){ dVector[1] }, "pair "
            "motif distance");
      }
      else {

        outputFile.printMetaLine((std::vector<int>){ window }, "window");
        outputFile.printMetaLine(motifPositions[0], "pair motif");
        outputFile.printMetaLine(dVector, "distance");
      }

      outputFile.close();
    }
    catch (int e) {

      if (e == EXIT_FAILURE)
        exit(e);
      else {

        std::cerr << "ERROR: Something unexpected happened!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    catch (...) {

      std::cerr << "ERROR: Something unexpected happened!" << std::endl;
      exit(EXIT_FAILURE);
    }

    running = false;
  }
}

void TsgGui::plotMotif() {

  if (!running) {

    running = true;
    double min, max;

    if (motifList.currentRow() != selMotif) {

      selMotif = motifList.currentRow();


      //upate subsequence marker
      upperSeries.clear();
      lowerSeries.clear();

      upperSeries.append(motifPositions[0][selMotif],
          100000000.0);
      upperSeries.append(motifPositions[0][selMotif] + window - 1,
          100000000.0);
      lowerSeries.append(motifPositions[0][selMotif],
          -100000000.0);
      lowerSeries.append(motifPositions[0][selMotif] + window - 1,
          -100000000.0);


      //update motif location plot
      motifSeries.clear();

      motifSeries.append(motifPositions[0][selMotif],
          timeSeries[motifPositions[0][selMotif]]);
      min = timeSeries[motifPositions[0][selMotif]];
      max = min;

      for (int i = motifPositions[0][selMotif] + 1;
          i < motifPositions[0][selMotif] + window; i++) {

        if (timeSeries[i] < min)
          min = timeSeries[i];
        if (timeSeries[i] > max)
          max = timeSeries[i];
        motifSeries.append(i, timeSeries[i]);
      }

      motifChart.axes(Qt::Horizontal)[0]->setRange(motifPositions[0][selMotif],
          motifPositions[0][selMotif] + window - 1);
      motifChart.axes(Qt::Vertical)[0]->setRange(min - 0.1 * (max - min), max
          + 0.1 * (max - min));
    }

    running = false;
  }
}

void TsgGui::showGenHelp() {

  genMes.show();
}

void TsgGui::showTypeHelp() {

  typeMes.show();
}

void TsgGui::showMethodHelp() {

  methodMes.show();
}

void TsgGui::showLengthHelp() {

  lengthMes.show();
}

void TsgGui::showWindowHelp() {

  windowMes.show();
}

void TsgGui::showSizeHelp() {

  sizeMes.show();
}

void TsgGui::showNoiseHelp() {

  noiseMes.show();
}

void TsgGui::showDeltaHelp() {

  deltaMes.show();
}

void TsgGui::showHeightHelp() {

  heightMes.show();
}

void TsgGui::showStepHelp() {

  stepMes.show();
}

void TsgGui::showTimesHelp() {

  timesMes.show();
}

void TsgGui::showMaxiHelp() {

  maxiMes.show();
}

