///\file tsm.cpp
///
///\brief File contains a top set motif discovery algorithm.
///
///This algorithm computes the top set motif of a time series by evaluating
///iteratively the diagonals of the distance matrix of the time series.

#include <tsm.hpp>


namespace tsg {

  TSM::TSM(const rseq &timeSeries_in, const rseq &sums_in, const rseq
      &sumSquares_in) : timeSeries(timeSeries_in), sums(sums_in),
      sumSquares(sumSquares_in) {

    if (timeSeries_in.empty()) {

      std::cerr << "ERROR: Time series is empty!" << std::endl;
      throw(EXIT_FAILURE);
    }

    if (sums_in.empty()) {

      std::cerr << "ERROR: Running sums are empty!" << std::endl;
      throw(EXIT_FAILURE);
    }

    if (sumSquares_in.empty()) {

      std::cerr << "ERROR: Running sum of squares are empty!" << std::endl;
      throw(EXIT_FAILURE);
    }
  }

  void TSM::zNormalPAA(rseq &paa_out, const int pos_in, const int window_in,
      const int paa_in, const bool roll_in) {

    if (!paa_out.empty()) {

      paa_out.clear();
      paa_out.resize(0);
    }

    double rW = 1.0 / (double)window_in;
    double c = (double)window_in / (double)paa_in;
    double rC = (double)paa_in / (double)window_in;
    double s = 0.0;
    double mean = rW * sums[pos_in];
    double sig = (rW * (sumSquares[pos_in] - rW * (sums[pos_in]
            * sums[pos_in])));

    if (sig < 1.0)
      sig = 1.0;
    else
      sig = sqrt(sig);

    double rSig = 1.0 / sig;
    double rCrS = rC * rSig;

    for (int i = 0; i < paa_in; i++) {

      s = 0;

      int start = c * i;
      int end = (int)(c * i + c);

      //rolling case, use last word
      if (roll_in) {

        s = lastPAA[i] + rCrS * (timeSeries[pos_in - 1 + end]
            - timeSeries[pos_in - 1 + start]);
      }
      //regular case, incrementally from time series
      else {

        for (int j = start; j < end; j++) {

            s += (timeSeries[pos_in + j] - mean) * rSig;
            s *= rC;
        }
      }

      paa_out.push_back(s);
    }

    lastPAA = paa_out;
  }

  double TSM::invNormalCDF(const double p_in, const double prec_in) {

    if (p_in <= 0.0)
      return -std::numeric_limits<double>::infinity();

    if (p_in >= 1.0)
      return std::numeric_limits<double>::infinity();

    double p = 0.0;
    double err = 1.0;
    double conq = 0.25;
    double inv = 0.5;
    double sqrtH = sqrt(0.5);
    double last = -1.0;

    while (fabs(err) > prec_in) {

      p = 0.5 * erfc(-inv * sqrtH);
      err = p_in - p;

      if (p_in == p)
        return inv;
      else if (p_in < p) {

        inv -= conq;

        if (last > 0.0) {

          conq *= 0.5;
          last *= -1.0;
        }
      }
      else {

        inv += conq;

        if (last < 0.0) {

          conq *= 0.5;
          last *= -1.0;
        }
      }
    }

    return inv;
  }

  void TSM::invNormalCDF(const int alphabet_in, rseq &p_out, const double
      prec_in) {

    if (!p_out.empty()) {

      p_out.clear();
      p_out.resize(0);
    }

    double n = alphabet_in;
    double rN = 1.0 / n;

    for (int i = 0; i < (int)n + 1; i++)
      p_out.push_back(invNormalCDF(rN * i, prec_in));
  }

  void TSM::zNormalSAX(word &sax_out, const int pos_in, const int window_in,
      const int paa_in, const rseq &breakpoints_in, const bool roll_in) {

    if (!sax_out.empty()) {

      sax_out.clear();
      sax_out.resize(0);
    }

    rseq paa;

    zNormalPAA(paa, pos_in, window_in, paa_in, roll_in);

    for (int i = 0; i < (int)paa.size(); i++)
      for (int j = 0; j < (int)breakpoints_in.size(); j++)
        if (paa[i] < breakpoints_in[j]) {

          sax_out.push_back(j);
          j = (int)breakpoints_in.size();
        }
  }

  double TSM::saxDist(const word &word0_in, const word word1_in, const int
      window_in, const rseq &breakpoints_in) {

    double dist = std::numeric_limits<double>::infinity();

    if (word0_in.size() != word1_in.size() || word0_in.empty())
      return dist;

    double m = (double)word0_in.size();
    double saxDist = 0;

    dist = 0.0;

    for (int i = 0; i < m; i++) {

      if (abs(word0_in[i] - word1_in[i]) > 1) {

        if (word0_in[i] < word1_in[i])
          saxDist = breakpoints_in[word0_in[i]]
            - breakpoints_in[word1_in[i] - 1];
        else
          saxDist = breakpoints_in[word1_in[i]]
            - breakpoints_in[word0_in[i] - 1];

        dist += saxDist * saxDist;
      }
    }

    return sqrt((double)window_in * dist / m);
  }

  double TSM::dist(const int pos0_in, const int pos1_in, const int window_in) {

    int window = window_in;
    int length = (int)timeSeries.size();

    if (pos0_in + window > length || pos0_in < 0) {

      std::cerr << "ERROR: Position of first subsequence " <<
        pos0_in << " in similarity function is wrong!" << std::endl;
      throw(EXIT_FAILURE);
    }

    if (pos1_in + window > length || pos1_in < 0) {

      std::cerr << "ERROR: Position of second subsequence " <<
        pos1_in << " in similarity function is wrong!" << std::endl;
      throw(EXIT_FAILURE);
    }

    double rWindow = 1.0 / (double)window;

    double mean0 = sums[pos0_in] * rWindow;
    double stdDev0 = sumSquares[pos0_in] * rWindow  - mean0 * mean0;
    stdDev0 = stdDev0 < 1.0 ? 1.0 : sqrt(stdDev0);

    double mean1 = sums[pos1_in] * rWindow;
    double stdDev1 = sumSquares[pos1_in] * rWindow  - mean1 * mean1;
    stdDev1 = stdDev1 < 1.0 ? 1.0 : sqrt(stdDev1);

    //compute the distance
    double dist = 0.0;
    double norm0;
    double norm1;
    double diff;

    for (int i = 0; i < window; i++) {

      norm0 = (timeSeries[pos0_in + i] - mean0) / stdDev0;
      norm1 = (timeSeries[pos1_in + i] - mean1) / stdDev1;
      diff = norm0 - norm1;
      dist += diff * diff;
    }

    return sqrt(dist);
  }

  double TSM::distADM(const int i_in, const int j_in) {

    int i = i_in;
    int j = j_in;

    if (i == j)
      return 0.0;

    if (j > i) {

      i = j_in;
      j = i_in;
    }

    return admDist[i - 1][j];
  }

  void TSM::adm(const iseq &neighborhood_in, const int window_in, const double
      range_in) {

    if (admDist.empty()) {

      admDist.resize(neighborhood_in.size() - 1);
    }
    else {

      for (auto &item : admDist) {

        item.clear();
        item.resize(0);
      }

      admDist.clear();
      admDist.resize(neighborhood_in.size() - 1);
    }

    if ((int)admDist.size() > 0) {

      //maximum distance matrix initialized to 0.0
      rseqs &a = admDist;
      //minimum distance matrix initialized to infinity
      rseqs m(a.size(), rseq(0));

      for (int i = 0; i < (int)a.size(); i++) {

        for (int j = 0; j <= i; j++) {

          a[i].push_back(0.0);
          m[i].push_back(std::numeric_limits<double>::infinity());
        }
      }

      std::mt19937 randomEngine(std::random_device().entropy()
          ? std::random_device()()
          : (unsigned int)
          (std::chrono::system_clock::now().time_since_epoch().count()));

      std::uniform_int_distribution<int> distributionX(0, (int)a.size() - 1);

      int x = 0;
      int y = 0;
      double d = 0.0;

      //compute exact distance for random number of cells
      for (int i = 0; i < 2 * (int)a.size(); i++) {

        x = distributionX(randomEngine);

        std::uniform_int_distribution<int> distributionY(0, x);

        y = distributionY(randomEngine);

        d = dist(neighborhood_in[x + 1], neighborhood_in[y], window_in);
        a[x][y] = d;
        m[x][y] = d;
      }

      for (int k = 0; k < (int)a.size(); k++) {

        for (int i = 0; i <= k; i++) {

          for (int j = 0; j <= i; j++) {

            a[i][j] = std::max(std::max(a[i][j], a[k][i] - m[k][j]), a[k][j]
                - m[k][i]);
            m[i][j] = std::min(m[i][j], m[k][i] + m[k][j]);
          }
        }
      }

      for (int i = 0; i < (int)a.size(); i++)
        for (int j = 0; j <= i; j++)
          if (a[i][j] < range_in)
            a[i][j] = dist(neighborhood_in[i + 1], neighborhood_in[j], window_in);
    }
  }

  int TSM::tsm(iseq &motif_out, const int window_in, const double range_in,
      const int paa_in, const rseq &breakpoints_in) {

    double range = range_in;

    if (range < 0.0)
      range = -range;

    if (!motif_out.empty()) {

      motif_out.clear();
      motif_out.resize(0);
    }

    int window = window_in;
    int length = (int)(timeSeries.size());

    //compute buckets containing all indices of subsequences with similar
    //z-normalized sax representation
    std::vector<std::pair<word, iseq>> buckets;
    std::vector<std::pair<word, iseq>>::iterator itB;
    word sax = "";

    zNormalSAX(sax, 0, window, paa_in, breakpoints_in);

    itB = find_if(buckets.begin(), buckets.end(), [&sax](const
          std::pair<word, iseq> &e){return e.first == sax;});

    if (itB == buckets.end())
      buckets.push_back(std::pair<word, iseq>(word(sax), { 0 }));
    else
      itB->second.push_back(0);

    for (int i = 1; i < length - window + 1; i++) {

      zNormalSAX(sax, i, window, paa_in, breakpoints_in, true);

      itB = find_if(buckets.begin(), buckets.end(), [&sax](const
            std::pair<word, iseq> &e){return e.first == sax;});

      if (itB == buckets.end())
        buckets.push_back(std::pair<word, iseq>(word(sax), { i }));
      else
        itB->second.push_back(i);
    }

    //sort in descending bucket size order
    std::sort(buckets.begin(), buckets.end(), [](std::pair<word, iseq> const &a,
          std::pair<word, iseq> const &b) {
        return a.second.size() > b.second.size();
        });

    //build neighborhood
    iset collect;
    iseq neighborhood;

    for (auto &item : buckets[0].second)
      collect.insert(item);

    for (int i = 1; i < (int)buckets.size(); i++)
      if (saxDist(buckets[0].first, buckets[i].first, window_in,
            breakpoints_in) < range)
        for (int j = 0; j < (int)buckets[i].second.size(); j++)
          collect.insert(buckets[i].second[j]);

    for (auto &item : collect)
      neighborhood.push_back(item);

    std::sort(neighborhood.begin(), neighborhood.end());

    collect.clear();

    //most promissing candidate
    iseq &mpc = motif_out;
    iseq motif;

    //select the largest unprocessed bucket
    for (int i = 0; i < (int)buckets.size() &&
        neighborhood.size() > mpc.size(); i++) {

      //compute the distance matrix of the neighborhood
      adm(neighborhood, window_in, range);

      //get largest set motif in the neighborhood
      for (int j = 0; j < (int)neighborhood.size(); j++) {

        //if there does not exist a set motif in the remaining neighborhood
        //matching more subsequences then the most promissing one
        if ((int)neighborhood.size() - j <= (int)mpc.size())
          break;

        if (!motif.empty()) {

          motif.clear();
          motif.resize(0);
        }

        //add the motif subsequence
        motif.push_back(neighborhood[j]);

        //get all subsequence matching the motif subsequence
        for (int k = 0; k < (int)neighborhood.size(); k++)
          if (distADM(j, k) < range)
            motif.push_back(neighborhood[k]);
        //motif is sorted in ascending order as the neighborhood is already
        //sorted

        //filter overlapping subsequences from motif
        int last = 0;

        for (int k = 1; k < (int)motif.size(); k++) {

          if (abs(motif[last] - motif[k]) < window_in) {

            //delete overlapping subsequences
            motif.erase(motif.begin() + k);
            k--;
          }
          else
            last = k;
        }

        //check if motif is better
        if (mpc.size() < motif.size()) {

          mpc.clear();

          for (int k = 0; k < (int)motif.size(); k++)
            mpc.push_back(motif[k]);
        }
      }

      //build next neighborhood
      neighborhood.clear();
      neighborhood.resize(0);

      for (auto &item : buckets[i].second)
        collect.insert(item);

      for (int j = 0; j < (int)buckets.size(); j++)
        if (saxDist(buckets[j].first, buckets[j].first, window_in,
              breakpoints_in) < range)
          for (int k = 0; k < (int)buckets[j].second.size(); k++)
            collect.insert(buckets[j].second[k]);

      for (auto &item : collect)
        neighborhood.push_back(item);

      std::sort(neighborhood.begin(), neighborhood.end());

      collect.clear();
    }

    return (int)mpc.size();
  }
}
