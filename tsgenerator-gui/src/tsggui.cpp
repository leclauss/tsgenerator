///\file tsggui.cpp
///
///\brief Source file of the TSGenerator GUI.
///
///This is the source file of the TSGenerator GUI.

#include <tsggui.hpp>


TsgGui::TsgGui(int &argc, char *argv[]) : QApplication(argc, argv),
  upperSeriess(1000), lowerSeriess(1000), markers(1000), motifSeries(1000) {

  //enable to select multiple lines in list
  motifList.setSelectionMode(QAbstractItemView::ExtendedSelection);

  //set initialize widgets
  fileMenu.setTitle("File");
  openAct.setText("open");
  openFileDialog.setWindowTitle("open file");
  openFileDialog.setFileMode(QFileDialog::Directory);
  openFileDialog.setOption(QFileDialog::ShowDirsOnly, true);
  openFileDialog.setDirectory(".");
  openFileDialog.setNameFilter("*.csv");
  saveAct.setText("save");
  infoMenu.setTitle("Info");
  copywriteAct.setText("license");
  copywriteMes.setWindowTitle("license");
  copywriteMes.setText("MIT License\n\n"
      "Copyright (c) 2018 Rafael Moczalla\n\n"
      "Permission is hereby granted, free of charge, to any person obtaining "
      "a copy of this software and associated documentation files (the "
      "\"Software\"), to deal in the Software without restriction, including "
      "without limitation the rights to use, copy, modify, merge, publish, "
      "distribute, sublicense, and/or sell copies of the Software, and to "
      "permit persons to whom the Software is furnished to do so, subject to "
      "the following conditions:\n\n"
      "The above copyright notice and this permission notice shall be "
      "included in all copies or substantial portions of the Software.\n\n"
      "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, "
      "EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF "
      "MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. "
      "IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY "
      "CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, "
      "TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE "
      "SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.");
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
  startButton.setText("start");
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
  distLabel.setText("distance");
  rangeLabel.setText("range");
  browseLabel.setText("data file:");
  browseButton.setText("browse");

  //setup menu actions
  connect(&openAct, SIGNAL(triggered()), this, SLOT(openData()));
  connect(&saveAct, SIGNAL(triggered()), this, SLOT(saveData()));

  connect(&copywriteAct, SIGNAL(triggered()), this, SLOT(showCopywrite()));

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
  fileMenu.addAction(&openAct);
  fileMenu.addAction(&saveAct);
  menuBar.addMenu(&fileMenu);

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

  infoMenu.addAction(&copywriteAct);
  menuBar.addMenu(&infoMenu);

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

  //setup save button action
  connect(&startButton, SIGNAL(clicked()), this, SLOT(generateTS()));
  connect(&browseButton, SIGNAL(clicked()), this, SLOT(openData()));

  //set subsqeuence marker
  motifPen.setColor(0xffbb00);
  brush.setColor(0xffddaa);
  brush.setStyle(Qt::SolidPattern);

  for (int i = 0; i < (int)lowerSeriess.size(); i++) {

    markers[i].setPen(motifPen);
    markers[i].setBrush(brush);
    markers[i].setUpperSeries(&upperSeriess[i]);
    markers[i].setLowerSeries(&lowerSeriess[i]);
  }

  //set time series color
  tsPen.setColor(0x0076ff);
  tsSeries.setPen(tsPen);

  //set motif plot color
  for (int i = 0; i < (int)motifSeries.size(); i++)
    motifSeries[i].setPen(motifPen);

  //add time series plot
  for (int i = 0; i < (int)lowerSeriess.size(); i++)
    tsChart.addSeries(&markers[i]);

  tsChart.addSeries(&tsSeries);
  tsChart.createDefaultAxes();
  tsChart.legend()->hide();

  tsChartView.setChart(&tsChart);
  tsChartView.setRenderHint(QPainter::Antialiasing);
  tsChartView.setMinimumSize(1200, 200);
  tsChartView.setMaximumSize(30000, 500);

  //add motif location plot
  for (int i = 0; i < (int)motifSeries.size(); i++)
    motifChart.addSeries(&motifSeries[i]);

  motifChart.createDefaultAxes();
  motifChart.legend()->hide();

  motifChartView.setChart(&motifChart);
  motifChartView.setRenderHint(QPainter::Antialiasing);
  motifChartView.setMinimumSize(20, 250);
  motifChartView.setMaximumSize(30000, 600);

  //add plot motif location trigger
  connect(&motifList, SIGNAL(itemSelectionChanged()), this, SLOT(plotMotif()));

  //add open browser when enter is pressed
  connect(&browseText, SIGNAL(returnPressed()), this, SLOT(openData()));

  //layout all items
  browseLayout.addWidget(&browseLabel);
  browseLayout.addWidget(&browseText);
  browseLayout.addWidget(&browseButton);

  browseBox.setTitle("data source");
  browseBox.setLayout(&browseLayout);

  tsLayout.addWidget(&tsChartView);

  tsBox.setTitle("time series");
  tsBox.setLayout(&tsLayout);

  motifLayout.addWidget(&motifChartView);

  motifBox.setTitle("z-normalized zoom");
  motifBox.setLayout(&motifLayout);

  motifListLayout.addWidget(&motifList);

  motifListBox.setTitle("motif locations");
  motifListBox.setLayout(&motifListLayout);

  layout.addWidget(&browseBox, 0, 0, 1, 6);
  layout.addWidget(&genLabel, 1, 0);
  layout.addWidget(&genDrop, 2, 0);
  layout.addWidget(&typeLabel, 1, 1);
  layout.addWidget(&typeDrop, 2, 1);
  layout.addWidget(&methodLabel, 1, 2);
  layout.addWidget(&methodDrop, 2, 2);
  layout.addWidget(&lengthLabel, 1, 3);
  layout.addWidget(&lengthText, 2, 3);
  layout.addWidget(&windowLabel, 1, 4);
  layout.addWidget(&windowText, 2, 4);
  layout.addWidget(&sizeLabel, 1, 5);
  layout.addWidget(&sizeText, 2, 5);
  layout.addWidget(&noiseLabel, 3, 0);
  layout.addWidget(&noiseText, 4, 0);
  layout.addWidget(&deltaLabel, 3, 1);
  layout.addWidget(&deltaText, 4, 1);
  layout.addWidget(&heightLabel, 3, 2);
  layout.addWidget(&heightText, 4, 2);
  layout.addWidget(&stepLabel, 3, 3);
  layout.addWidget(&stepText, 4, 3);
  layout.addWidget(&timesLabel, 3, 4);
  layout.addWidget(&timesText, 4, 4);
  layout.addWidget(&maxiLabel, 3, 5);
  layout.addWidget(&maxiText, 4, 5);
  layout.addWidget(&startButton, 5, 5);
  layout.addWidget(&tsBox, 6, 0, 1, -1);
  layout.addWidget(&motifBox, 7, 0, 8, 3);
  layout.addWidget(&motifListBox, 7, 3, 8, 2);
  layout.addWidget(&distLabel, 7, 5);
  layout.addWidget(&distText, 8, 5);
  layout.addWidget(&rangeLabel, 9, 5);
  layout.addWidget(&rangeText, 10, 5);

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

void TsgGui::loadData() {

  if (!dVector.empty()) {

    //update smallest distance and motif range
    if (gen == "pair motif") {

      distText.setText(std::to_string(dVector[0]).c_str());
      rangeText.setText("");
    }
    else {

      distText.setText(std::to_string(dVector[1]).c_str());
      rangeText.setText(std::to_string(dVector[0]).c_str());
    }
  }

  //print motif select list
  motifLocs.clear();
  motifList.clear();

  int start = 0;

  if (gen == "set motif") {

    if (!motifPositions.empty() && !motifPositions[0].empty()) {

      std::stringstream setMotif;
      setMotif << "set motif loc: " << motifPositions[0][0];
      motifLocs.append(setMotif.str().c_str());

      start = 1;
    }
  }
  else if (gen == "latent motif") {

    if (!motif.empty() && !motif[0].empty()) {

      std::stringstream latent;
      latent << "latent motif";
      motifLocs.append(latent.str().c_str());
    }
  }
  else {

    if (!motifPositions.empty() && !motifPositions[0].empty()) {

      //first pair motif sequence
      std::stringstream pair0;
      pair0 << "pair motif loc: " << motifPositions[0][0];
      motifLocs.append(pair0.str().c_str());

      //second pair motif sequence
      std::stringstream pair1;
      pair1 << "pair motif loc: " << motifPositions[0][1];
      motifLocs.append(pair1.str().c_str());

      start = 2;
    }
  }

  if (!motifPositions.empty() && !motifPositions[0].empty()) {

    for (int i = start; i < (int)motifPositions[0].size(); i++) {

      std::stringstream idx;
      idx << "matching loc: " << motifPositions[0][i];
      motifLocs.append(idx.str().c_str());
    }
  }

  //add algorithm results to motif list
  tsg::word alg;

  for (auto &item : disc) {

    alg = item.first;

    for (int i = 0; i < (int)item.second.size(); i++)
      motifLocs.append((alg + " > matchings " + std::to_string(i)).c_str());
  }

  //print the list
  motifList.addItems(motifLocs);
  motifList.setCurrentRow(0);

  //highlight selected subsequences and plot zooms
  QList selected = motifList.selectedItems();
  int j = 0;
  int pos;

  double min = std::numeric_limits<double>::infinity();
  double max = -std::numeric_limits<double>::infinity();

  for (int i = 0; i < (int)markers.size(); i++)
    markers[i].hide();

  for (int i = 0; i < (int)motifSeries.size(); i++)
    motifSeries[i].clear();

  double rWindow = 1.0 / (double)window;

  for (int i = 0; i < selected.size(); i++) {

    if (j < (int)lowerSeriess.size()) {

      double mean = 0.0;
      double std = 0.0;
      double value = 0.0;

      if (selected[i]->text().indexOf(": ") >= 0)
        {

        pos = selected[i]->text().split(": ")[1].toInt();

        upperSeriess[j].clear();
        lowerSeriess[j].clear();

        upperSeriess[j].append(pos, 100000000.0);
        upperSeriess[j].append(pos + window - 1, 100000000.0);
        lowerSeriess[j].append(pos, -100000000.0);
        lowerSeriess[j].append(pos + window - 1, -100000000.0);

        markers[j].show();

        if (!timeSeries.empty() && !motifPositions.empty() &&
            !motifPositions[0].empty()) {

          for (int i = pos; i < pos + window; i++)
            mean += timeSeries[i];

          mean *= rWindow;

          for (int i = pos; i < pos + window; i++)
            std += timeSeries[i] * timeSeries[i] - mean * mean;

          std *= rWindow;

          if (std <= 1.0)
            std = 1.0;
          else
            std = sqrt(std);

          std = 1.0 / std;

          for (int i = 0; i < window; i++) {

            value = (timeSeries[pos + i] - mean) * std;

            if (value < min)
              min = value;
            if (value > max)
              max = value;

            motifSeries[j].append(i, value);
          }

          j++;
        }
      }
      else if (QString::compare(selected[i]->text(), "latent motif",
            Qt::CaseInsensitive) == 0) {

        if (!motif.empty() && !motif[0].empty()) {

          for (int i = 0; i < window; i++)
            mean += motif[0][i];

          mean *= rWindow;

          for (int i = 0; i < window; i++)
            std += motif[0][i] * motif[0][i] - mean * mean;

          std *= rWindow;

          if (std <= 1.0)
            std = 1.0;
          else
            std = sqrt(std);

          std = 1.0 / std;

          for (int k = 0; k < (int)motif[0].size(); k++) {

            value = (motif[0][k] - mean) * std;

            if (value < min)
              min = value;
            if (value > max)
              max = value;

            motifSeries[j].append(k, value);
          }

          j++;
        }
      }
    }
  }

  if (selected.size() > 0 && !motifPositions.empty() &&
    !motifPositions[0].empty()) {

    motifChart.axes(Qt::Horizontal)[0]->setRange(0, window - 1);
    motifChart.axes(Qt::Vertical)[0]->setRange(min - 0.1 * (max - min), max
        + 0.1 * (max - min));
  }

  //plot ts
  tsSeries.clear();

  if (!timeSeries.empty()) {

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

    tsChart.axes(Qt::Vertical)[0]->setRange(min - 0.1 * (max - min), max
        + 0.1 * (max - min));

    tsChart.axes(Qt::Horizontal)[0]->setRange(0.0, length);
  }
}

void TsgGui::generateTS() {

  if (!running) {

    running = true;
    bool success = false;
    tsg::word in;

    //read the input
    in = lengthText.text().toStdString();
    if (in.length() > 0 && in.find_first_not_of("0123456789") == in.npos)
      length = std::stoi(in);

    in = windowText.text().toStdString();
    if (in.length() > 0 && in.find_first_not_of("0123456789") == in.npos)
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
    if (in.length() > 0 && in.find_first_not_of("0123456789") == in.npos)
      motifSize = std::stoi(in);

    in = heightText.text().toStdString();
    if (check_if_float(in))
      height = std::stod(in);

    in = stepText.text().toStdString();
    if (in.length() > 0 && in.find_first_not_of("0123456789") == in.npos)
      step = std::stoi(in);

    in = timesText.text().toStdString();
    if (in.length() > 0 && in.find_first_not_of("0123456789") == in.npos)
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
    catch (...) {

      std::cerr << "ERROR: Something unexpected happened!" << std::endl;
    }

    if (success) {

      //clear discovered motifs data
      disc.clear();

      loadData();
    }
    else {
      std::cerr << "ERROR: Couldn't create time series." << std::endl;
    }

    running = false;
  }
}

void TsgGui::openData() {

  if (!running) {

    running = true;
    size_t pos;

    QString path = browseText.text();
    tsg::word tmp = path.toStdString();

    if (!tmp.empty()) {

      pos = tmp.find_last_of("/");

      if (pos == tmp.npos) {

        openFileDialog.setDirectory("/");
        openFileDialog.selectFile(tmp.c_str());
      }
      else {

        openFileDialog.setDirectory(tmp.substr(0, pos).c_str());
        openFileDialog.selectFile(tmp.substr(pos + 1).c_str());
      }
    }

    openFileDialog.exec();
    QStringList selectedFiles = openFileDialog.selectedFiles();

    if (selectedFiles.isEmpty() ||
        !std::filesystem::exists(selectedFiles[0].toStdString())) {

      running = false;
      return;
    }

    tsg::word dirPath = selectedFiles[0].toStdString();
    tsg::par outPaths;
    tsg::word tsFilePath;
    tsg::word metaFilePath;
    tsg::word posIntChars("0123456789");

    //process all files in the directory
    for (const auto &entry : std::filesystem::directory_iterator(dirPath)) {

      //clear discovered motifs data
      disc.clear();

      tmp = entry.path().filename();

      pos = tmp.rfind(".out");

      //motif discovery output files
      if (pos != tmp.npos) {

        outPaths.push_back(entry.path().filename());
        continue;
      }

      pos = tmp.rfind(".csv");

      if (pos != tmp.npos) {

        tmp.replace(pos, 4, "");
        pos = tmp.rfind("time_series_meta_");

        //meta data file
        if (pos != tmp.npos) {

          tmp.replace(pos, 17, "");

          if (tmp.find_first_not_of(posIntChars) == tmp.npos)
            metaFilePath = dirPath + "/time_series_meta_" + tmp + ".csv";

          continue;
        }

        pos = tmp.rfind("time_series_");

        //time series file
        if (pos != tmp.npos) {

          tmp.replace(pos, 12, "");

          if (tmp.find_first_not_of(posIntChars) == tmp.npos)
            tsFilePath = dirPath + "/time_series_" + tmp + ".csv";
        }
      }
    }

    browseText.setText(dirPath.c_str());

    std::ifstream tsFile;
    std::ifstream metaFile;
    double value;
    tsg::word doubleChars("-.0123456789");
    tsg::word intChars("-0123456789");
    tsg::word alphabetChars(" abcdefghijklmnopqrstuvwxyz"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ");

    if (std::filesystem::exists(tsFilePath)) {

      if (!timeSeries.empty())
        timeSeries.clear();

      tsFile.open(tsFilePath);

      while (doubleChars.find(tsFile.peek()) == doubleChars.npos &&
          tsFile.peek() != EOF)
        tsFile.ignore();

      while (tsFile.peek() != EOF && tsFile >> value) {

        timeSeries.push_back(value);

        while (doubleChars.find(tsFile.peek()) == doubleChars.npos &&
            tsFile.peek() != EOF)
          tsFile.ignore();
      }

      length = (int)timeSeries.size();

      tsFile.close();
    }

    if (std::filesystem::exists(tsFilePath) &&
        std::filesystem::exists(metaFilePath)) {

      if (!motif.empty()) {

        for (auto &item : motif)
          if (!item.empty())
            item.clear();

        motif.clear();
      }

      motif.resize(1);

      if (!motifPositions.empty()) {

        for (auto &item : motifPositions)
          if (!item.empty())
            item.clear();

        motifPositions.clear();
      }

      motifPositions.resize(2);

      if (!dVector.empty()) {

        dVector.clear();
      }

      metaFile.open(metaFilePath);

      //read window
      while (posIntChars.find(metaFile.peek()) == posIntChars.npos &&
          metaFile.peek() != EOF)
        metaFile.ignore();

      if (metaFile.peek() != EOF)
        metaFile >> window;

      while (alphabetChars.find(metaFile.peek()) == alphabetChars.npos &&
          metaFile.peek() != EOF)
        metaFile.ignore();

    std::cout << "file: " << metaFilePath << std::endl;
      //check which motif generator
      if (metaFile.peek() == 's') {

        gen = "set motif";

        while (posIntChars.find(metaFile.peek()) == posIntChars.npos &&
            metaFile.peek() != EOF)
          metaFile.ignore();

        while (alphabetChars.find(metaFile.peek()) == alphabetChars.npos &&
            metaFile.peek() != EOF)
          metaFile.ignore();
      }
      else if (metaFile.peek() == 'l') {

        gen = "latent motif";

        while (doubleChars.find(metaFile.peek()) == doubleChars.npos &&
            metaFile.peek() != EOF)
          metaFile.ignore();

        while (metaFile.peek() != '\n' && metaFile.peek() != EOF &&
            metaFile >> value) {

          motif[0].push_back(value);

          while (doubleChars.find(metaFile.peek()) == doubleChars.npos &&
              metaFile.peek() != '\n' && metaFile.peek() != EOF)
            metaFile.ignore();
        }

        while (alphabetChars.find(metaFile.peek()) == alphabetChars.npos &&
            metaFile.peek() != EOF)
          metaFile.ignore();
      }
      else {

        gen = "pair motif";
      }

      //check if we have matching subsequences
      int pos;
      int start = 0;

      if (metaFile.peek() == 'm') {

        while (posIntChars.find(metaFile.peek()) == posIntChars.npos &&
            metaFile.peek() != EOF)
          metaFile.ignore();

        if (metaFile.peek() != EOF)
          metaFile >> pos;

        motifPositions[start].push_back(pos);

        while (posIntChars.find(metaFile.peek()) == posIntChars.npos &&
            metaFile.peek() != '\n' && metaFile.peek() != EOF)
          metaFile.ignore();

        //check if we have already a motif
        if (motif[0].empty())
          for (int i = 0; i < window; i++)
            motif[0].push_back(timeSeries[pos + i]);

        //get positions
        while (metaFile.peek() != '\n' && metaFile.peek() != EOF &&
            metaFile >> pos) {

          motifPositions[start].push_back(pos);

          if (pos + window - 1 > length)
            length = pos + window - 1;

          while (posIntChars.find(metaFile.peek()) == posIntChars.npos &&
              metaFile.peek() != '\n' && metaFile.peek() != EOF)
            metaFile.ignore();
        }

        start++;

        while (alphabetChars.find(metaFile.peek()) == alphabetChars.npos &&
            metaFile.peek() != EOF)
          metaFile.ignore();

        //get range
        while (doubleChars.find(metaFile.peek()) == doubleChars.npos &&
            metaFile.peek() != EOF)
          metaFile.ignore();

        if (metaFile.peek() != EOF)
          metaFile >> value;

        dVector.push_back(value);

        while (alphabetChars.find(metaFile.peek()) == alphabetChars.npos &&
            metaFile.peek() != EOF)
          metaFile.ignore();
      }

      //pair motif locations
      while (posIntChars.find(metaFile.peek()) == posIntChars.npos &&
          metaFile.peek() != EOF)
        metaFile.ignore();

      while (metaFile.peek() != '\n' && metaFile >> pos) {

        motifPositions[start].push_back(pos);

        if (pos + window - 1 > length)
          length = pos + window - 1;

        while (posIntChars.find(metaFile.peek()) == posIntChars.npos &&
            metaFile.peek() != '\n' && metaFile.peek() != EOF)
          metaFile.ignore();
      }

      while (alphabetChars.find(metaFile.peek()) == alphabetChars.npos &&
          metaFile.peek() != EOF)
        metaFile.ignore();

      //pair motif distance
      while (doubleChars.find(metaFile.peek()) == doubleChars.npos &&
          metaFile.peek() != EOF)
        metaFile.ignore();

      if (metaFile.peek() != EOF)
        metaFile >> value;

      dVector.push_back(value);

      metaFile.close();
    }

    //add motif discovery results
    std::ifstream outFile;
    tsg::word line;
    int i = 0;
    size_t mPos = 0;

    for (auto &file : outPaths) {

      tmp = dirPath + "/" + file;

      i = 0;
      tsg::iseqs motifs(1);

      //read file and add motifs
      outFile.open(tmp);

      while (std::getline(outFile, line)) {

        mPos = line.find_last_of("-");

        if (line.find_last_not_of(intChars) == line.npos &&
            (mPos == 0 || mPos == line.npos)) {

          motifs[i].push_back(std::stoi(line));
        }
        else {

          motifs.resize(motifs.size() + 1);
          i++;
        }
      }

      pos = tmp.find_last_of("/");

      if (pos != tmp.npos)
        tmp.replace(0, pos + 1, "");

      tmp.replace(tmp.rfind(".out"), 4, "");

      disc.insert(std::pair<tsg::word, tsg::iseqs>(tmp, motifs));

      outFile.close();
    }

    loadData();

    running = false;
  }
}

void TsgGui::saveData() {

  if (!running) {

    running = true;

    if (!timeSeries.empty() ||
        (0 < window && !motifPositions.empty() && !dVector.empty())) {

      try {

        //output stuff
        OutputGenerator outputFile;

        outputFile.open();
        if (!timeSeries.empty())
          outputFile.printTimeSeriesVertical(timeSeries);

        if (0 < window && !motifPositions.empty() && !dVector.empty()) {

          if (gen == "set motif") {

            outputFile.printMetaLine((tsg::iseq){ window }, "window");
            outputFile.printMetaLine((tsg::iseq){ motifPositions[0][0] },
                "set motif");
            outputFile.printMetaLine(motifPositions[0], "matchings");
            outputFile.printMetaLine((tsg::rseq){ dVector[0] }, "range");
            outputFile.printMetaLine(motifPositions[1], "pair motif");
            outputFile.printMetaLine((tsg::rseq){ dVector[1] }, "pair "
                "motif distance");
          }
          else if (gen == "latent motif") {

            outputFile.printMetaLine((tsg::iseq){ window }, "window");
            outputFile.printMetaLine(motif[0], "latent motif");
            outputFile.printMetaLine(motifPositions[0], "matchings");
            outputFile.printMetaLine((tsg::rseq){ dVector[0] }, "range");
            outputFile.printMetaLine(motifPositions[1], "pair motif");
            outputFile.printMetaLine((tsg::rseq){ dVector[1] }, "pair "
                "motif distance");
          }
          else {

            outputFile.printMetaLine((tsg::iseq){ window }, "window");
            outputFile.printMetaLine(motifPositions[0], "pair motif");
            outputFile.printMetaLine(dVector, "distance");
          }
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
    }

    running = false;
  }
}

void TsgGui::plotMotif() {

  if (!running) {

    running = true;
    double min = std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();

    //upate subsequence marker
    QList selected = motifList.selectedItems();
    int j = 0;
    int pos;

    for (int i = 0; i < (int)markers.size(); i++)
      markers[i].hide();

    for (int i = 0; i < (int)motifSeries.size(); i++)
      motifSeries[i].clear();

    double rWindow = 1.0 / (double)window;

    for (int i = 0; i < selected.size(); i++) {

      if (j < (int)lowerSeriess.size()) {

        double mean = 0.0;
        double std = 0.0;
        double value = 0.0;

        if (selected[i]->text().indexOf(": ") >= 0 &&
            j < (int)lowerSeriess.size()) {

          pos = selected[i]->text().split(": ")[1].toInt();

          upperSeriess[j].clear();
          lowerSeriess[j].clear();

          upperSeriess[j].append(pos, 100000000.0);
          upperSeriess[j].append(pos + window - 1, 100000000.0);
          lowerSeriess[j].append(pos, -100000000.0);
          lowerSeriess[j].append(pos + window - 1, -100000000.0);

          markers[j].show();

          if (!timeSeries.empty() && !motifPositions.empty() &&
              !motifPositions[0].empty()) {

            for (int i = pos; i < pos + window; i++)
              mean += timeSeries[i];

            mean *= rWindow;

            for (int i = pos; i < pos + window; i++)
              std += timeSeries[i] * timeSeries[i] - mean * mean;

            std *= rWindow;

            if (std <= 1.0)
              std = 1.0;
            else
              std = sqrt(std);

            std = 1.0 / std;

            for (int i = 0; i < window; i++) {

              value = (timeSeries[pos + i] - mean) * std;

              if (value < min)
                min = value;
              if (value > max)
                max = value;

              motifSeries[j].append(i, value);
            }

            j++;
          }
        }
        else if (selected[i]->text().indexOf(" > ") >= 0 &&
            j < (int)lowerSeriess.size()) {

          //determine the matching subsequences
          tsg::word alg = selected[i]->text().toStdString();
          pos = alg.find(" > matchings ");
          int k = std::stoi(alg.substr(pos + 13));

          alg.replace(pos, alg.size() - pos, "");

          tsg::iseq motif = disc[alg][k];

          //print the subsequences
          for (int &pos : motif) {

            upperSeriess[j].clear();
            lowerSeriess[j].clear();

            upperSeriess[j].append(pos, 100000000.0);
            upperSeriess[j].append(pos + window - 1, 100000000.0);
            lowerSeriess[j].append(pos, -100000000.0);
            lowerSeriess[j].append(pos + window - 1, -100000000.0);

            markers[j].show();

            if (!timeSeries.empty() && !motifPositions.empty() &&
                !motifPositions[0].empty()) {

              for (int i = pos; i < pos + window; i++)
                mean += timeSeries[i];

              mean *= rWindow;

              for (int i = pos; i < pos + window; i++)
                std += timeSeries[i] * timeSeries[i] - mean * mean;

              std *= rWindow;

              if (std <= 1.0)
                std = 1.0;
              else
                std = sqrt(std);

              std = 1.0 / std;

              for (int i = 0; i < window; i++) {

                value = (timeSeries[pos + i] - mean) * std;

                if (value < min)
                  min = value;
                if (value > max)
                  max = value;

                motifSeries[j].append(i, value);
              }

              j++;
            }
          }
        }
        else if (QString::compare(selected[i]->text(), "latent motif",
              Qt::CaseInsensitive) == 0) {

          if (!motif.empty() && !motif[0].empty()) {

            for (int i = 0; i < window; i++)
              mean += motif[0][i];

            mean *= rWindow;

            for (int i = 0; i < window; i++)
              std += motif[0][i] * motif[0][i] - mean * mean;

            std *= rWindow;

            if (std <= 1.0)
              std = 1.0;
            else
              std = sqrt(std);

            std = 1.0 / std;

            for (int k = 0; k < (int)motif[0].size(); k++) {

              value = (motif[0][k] - mean) * std;

              if (value < min)
                min = value;
              if (value > max)
                max = value;

              motifSeries[j].append(k, value);
            }

            j++;
          }
        }
      }
    }

    if (selected.size() > 0 && !motifPositions.empty() &&
        !motifPositions[0].empty()) {

      motifChart.axes(Qt::Horizontal)[0]->setRange(0, window - 1);
      motifChart.axes(Qt::Vertical)[0]->setRange(min - 0.1 * (max - min), max
          + 0.1 * (max - min));
    }


    running = false;
  }
}

void TsgGui::showCopywrite() {

  copywriteMes.show();
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

