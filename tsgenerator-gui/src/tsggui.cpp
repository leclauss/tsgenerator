///\file tsggui.cpp
///
///\brief Source file of the TSGenerator GUI.
///
///This is the source file of the TSGenerator GUI.

#include <tsggui.hpp>


TsgGui::TsgGui(int argc, char *argv[]) : QApplication(argc, argv) {

  //set initial widgets
  idxLabel.setText("Motif Locations");
  startButton.setText("start");
  saveButton.setText("save");
  lengthLabel.setText("length");
  lengthText.setPlaceholderText(to_string(length).c_str());
  windowLabel.setText("window");
  windowText.setPlaceholderText(to_string(window).c_str());
  sizeLabel.setText("size");
  sizeText.setPlaceholderText(to_string(motifSize).c_str());
  deltaLabel.setText("delta");
  deltaText.setPlaceholderText(to_string(delta).c_str());
  noiseLabel.setText("noise");
  noiseText.setPlaceholderText(to_string(noise).c_str());
  heightLabel.setText("height");
  heightText.setPlaceholderText(to_string(height).c_str());
  stepLabel.setText("step");
  stepText.setPlaceholderText(to_string(step).c_str());
  timesLabel.setText("times");
  timesText.setPlaceholderText(to_string(times).c_str());
  maxiLabel.setText("maxi");
  maxiText.setPlaceholderText(to_string(maxi).c_str());
  typeLabel.setText("type");
  methodLabel.setText("method");

  //add the types
  types.append("box");
  types.append("triangle");
  types.append("semicircle");
  types.append("trapezoid");
  types.append("positive flank");
  types.append("negative flank");
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

  //setup buttons
  connect(&startButton, SIGNAL(clicked()), this, SLOT(generateTS()));
  connect(&saveButton, SIGNAL(clicked()), this, SLOT(saveTS()));

  //add time series plot
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
  layout.addWidget(&typeLabel, 0, 0);
  layout.addWidget(&typeDrop, 1, 0);
  layout.addWidget(&methodLabel, 0, 1);
  layout.addWidget(&methodDrop, 1, 1);
  layout.addWidget(&lengthLabel, 0, 2);
  layout.addWidget(&lengthText, 1, 2);
  layout.addWidget(&windowLabel, 0, 3);
  layout.addWidget(&windowText, 1, 3);
  layout.addWidget(&sizeLabel, 0, 4);
  layout.addWidget(&sizeText, 1, 4);
  layout.addWidget(&noiseLabel, 0, 5);
  layout.addWidget(&noiseText, 1, 5);
  layout.addWidget(&deltaLabel, 2, 0);
  layout.addWidget(&deltaText, 3, 0);
  layout.addWidget(&heightLabel, 2, 1);
  layout.addWidget(&heightText, 3, 1);
  layout.addWidget(&stepLabel, 2, 2);
  layout.addWidget(&stepText, 3, 2);
  layout.addWidget(&timesLabel, 2, 3);
  layout.addWidget(&timesText, 3, 3);
  layout.addWidget(&maxiLabel, 2, 4);
  layout.addWidget(&maxiText, 3, 4);
  layout.addWidget(&startButton, 4, 4);
  layout.addWidget(&saveButton, 4, 5);
  layout.addWidget(&tsChartView, 5, 0, 1, -1);
  layout.addWidget(&motifChartView, 6, 0, 2, 3);
  layout.addWidget(&idxLabel, 6, 3, 1, 3);
  layout.addWidget(&motifList, 7, 3, 1, 3);

  pane.setLayout(&layout);

  //build window
  gui.setCentralWidget(&pane);
  gui.resize(800, 300);
  gui.setWindowTitle("tsgenerator");
  gui.show();
}

TsgGui::~TsgGui() { }

void TsgGui::generateTS() {

  if (!running) {

    running = true;
    bool first = false;
    bool success = false;
    double min, max;
    string in;

    //read the input
    in = lengthText.text().toStdString();
    if (in.length() > 0 && in.find_first_not_of("0123456789") ==
        string::npos)
      length = stoi(in);

    in = windowText.text().toStdString();
    if (in.length() > 0 && in.find_first_not_of("0123456789") ==
        string::npos)
      window = stoi(in);

    in = deltaText.text().toStdString();
    if (check_if_float(in))
      delta = stod(in);

    in = noiseText.text().toStdString();
    if (check_if_float(in))
      noise = stod(in);

    in = typeDrop.currentText().toStdString();
    if (in.length())
      type = in;

    in = sizeText.text().toStdString();
    if (in.length() > 0 && in.find_first_not_of("0123456789") ==
        string::npos)
      motifSize = stoi(in);

    in = heightText.text().toStdString();
    if (check_if_float(in))
      height = stod(in);

    in = stepText.text().toStdString();
    if (in.length() > 0 && in.find_first_not_of("0123456789") ==
        string::npos)
      step = stoi(in);

    in = timesText.text().toStdString();
    if (in.length() > 0 && in.find_first_not_of("0123456789") ==
        string::npos)
      times = stoi(in);

    in = methodDrop.currentText().toStdString();
    if (in.length())
      method = in;

    in = maxiText.text().toStdString();
    if (check_if_float(in))
      maxi = stod(in);

    try {

      //generate the time series
      tsg::TSGenerator tSGenerator(length, window, delta, noise, type,
          motifSize, height, step, times, method, maxi);
      tSGenerator.run(timeSeries, dVector, windows, motifPositions);

      success = true;
    }
    catch (int e) {

      if (e != EXIT_FAILURE)
        cerr << "ERROR: Something unexpected happened!" << endl;
    }
    catch (...) {

      cerr << "ERROR: Something unexpected happened!" << endl;
    }

    if (success) {

      //check if first time series
      if (motifs.empty())
        first = true;

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

      outputFile.printTimeSeriesVertical(timeSeries, dVector, windows,
          motifPositions);
    }
    catch (int e) {

      if (e == EXIT_FAILURE)
        exit(e);
      else {

        cerr << "ERROR: Something unexpected happened!" << endl;
        exit(EXIT_FAILURE);
      }
    }
    catch (...) {

      cerr << "ERROR: Something unexpected happened!" << endl;
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

