///\file scrimpplusplus.cpp
///
///\brief File contains the Scrimp++ algorithm.
///
///This is SCRIMP++.
///
///Details of the SCRIMP++ algorithm can be found at:
///(author information ommited for ICDM review),
///"SCRIMP++: Motif Discovery at Interactive Speeds", submitted to ICDM 2018.
///
///edited by Rafael Moczalla
//

#include <scrimpplusplus.hpp>


void scrimpPP(const vector<double> &timeSeries_in, int &pos0_out, int
    &pos1_out, int window_in, double stepSize_in)
{
  auto g = default_random_engine(random_device().entropy()
    ? random_device()()
    : (unsigned int)(chrono::system_clock::now().time_since_epoch().count()));

  int window = window_in;
  double rWindow = 1.0 / window;

  int stepSize = (int)(floor(stepSize_in*window));

  if (stepSize_in < 0.0)
    stepSize = (int)(floor(0.4*window));

  int timeSeriesLength = (int)(timeSeries_in.size());

  // set exclusion zone
  int exclusionZone = window / 4;

  // set Matrix Profile Length
  int ProfileLength = timeSeriesLength - window + 1;

  // preprocess, statistics, get the mean and standard deviation of every
  // subsequence in the time series
  vector<double> ASum;
  ASum.push_back(timeSeries_in[0]);
  for (int i = 1; i < window; i++)
    ASum[0] += timeSeries_in[i];
  for (int i = 0; i < timeSeriesLength - window; i++)
    ASum.push_back(ASum[i] - timeSeries_in[i] + timeSeries_in[window + i]);
  vector<double> ASumSq;
  ASumSq.push_back(timeSeries_in[0] * timeSeries_in[0]);
  for (int i = 1; i < window; i++)
    ASumSq[0] += timeSeries_in[i] * timeSeries_in[i];
  for (int i = 0; i < timeSeriesLength - window; i++)
    ASumSq.push_back(ASumSq[i] - timeSeries_in[i] * timeSeries_in[i]
        + timeSeries_in[i + window] * timeSeries_in[i + window]);
  vector<double> AMean;
  for (int i = 0; i < ProfileLength; i++)
    AMean.push_back(ASum[i] * rWindow);
  vector<double> ASigmaSq;
  for (int i = 0; i < ProfileLength; i++)
    ASigmaSq.push_back(ASumSq[i] * rWindow - AMean[i] * AMean[i]);
  vector<double> ASigma;
  for (int i = 0; i < ProfileLength; i++)
    ASigma.push_back(sqrt(ASigmaSq[i]));

  //Initialize Matrix Profile and Matrix Profile Index
  vector <double> profile;
  vector <int> profileIndex;
  for (int i=0; i<ProfileLength; i++)
  {
    profile.push_back(numeric_limits<double>::infinity());
    profileIndex.push_back(0);
  }

  //int fftsize = pow(2,ceil(log2(timeSeriesLength)));
  int fftsize = timeSeriesLength; //fftsize must be at least 2*window
  fftsize = fftsize > 2 * window ? fftsize : 2 * window;


  /*******************************PreSCRIMP***********************************/

  fftw_plan plan;
  fftw_complex* ATime = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)
      * fftsize);
  fftw_complex* AFreq = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)
      * fftsize);

  for (int i = 0; i < fftsize; i++)
  {
    ATime[i][1] = 0;
    if (i < timeSeriesLength)
      ATime[i][0] = timeSeries_in[i];
    else
      ATime[i][0] = 0;
  }

  plan = fftw_plan_dft_1d(fftsize, ATime, AFreq, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_free(ATime);
  fftw_destroy_plan(plan);

  fftw_complex* queryTime = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)
      * fftsize);
  fftw_complex* queryFreq = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)
      * fftsize);
  fftw_complex* AQueryTime = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)
      * fftsize);
  fftw_complex* AQueryFreq = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)
      * fftsize);

  //Sample subsequences with a fixed stepSize, then random shuffle their
  //computation order
  vector<int> idx;
  for (int i = 0; i < timeSeriesLength - window + 1; i += stepSize)
    idx.push_back(i);
  shuffle(idx.begin(), idx.end(), g);

  vector<double> query(window);

  for (int idx_i = 0; idx_i < (int) idx.size(); idx_i++)
  {
    int i = idx[idx_i];
    for (int j = 0; j < window; j++)
    {
      query[j] = timeSeries_in[i + j];
    }
    double queryMean = AMean[i];
    double queryStd = ASigma[i];

    for (int j = 0; j < fftsize; j++)
    {
      queryTime[j][1] = 0;

      if (j < window)
        queryTime[j][0] = query[window - j - 1];
      else
        queryTime[j][0] = 0;
    }

    plan = fftw_plan_dft_1d(fftsize, queryTime, queryFreq, FFTW_FORWARD,
        FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (int j = 0; j < fftsize; j++)
    {
      AQueryFreq[j][0] = AFreq[j][0] * queryFreq[j][0] - AFreq[j][1]
        * queryFreq[j][1];
      AQueryFreq[j][1] = AFreq[j][1] * queryFreq[j][0] + AFreq[j][0]
        * queryFreq[j][1];
    }

    plan = fftw_plan_dft_1d(fftsize, AQueryFreq, AQueryTime, FFTW_BACKWARD,
        FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    int exclusionZoneStart = i - exclusionZone;
    int exclusionZoneEnd = i + exclusionZone;
    double minimumDistance = numeric_limits<double>::infinity();
    int minimumDistanceIndex = 0;
    double distance;
    for (int j = 0; j < timeSeriesLength - window + 1; j++)
    {
      if ((j > exclusionZoneStart) && (j < exclusionZoneEnd))
        distance = numeric_limits<double>::infinity();
      else
      {
        distance = 2 * (window - (AQueryTime[window + j - 1][0]
              / fftsize - window * AMean[j] * queryMean) / (ASigma[j]
                * queryStd));
      }

      if (distance < minimumDistance)
      {
        minimumDistance = distance;
        minimumDistanceIndex = j;
      }

      if (distance < profile[j])
      {
        profile[j] = distance;
        profileIndex[j] = i;
      }
    }
    profile[i] = minimumDistance;
    profileIndex[i] = minimumDistanceIndex;

    int j = profileIndex[i];
    double lastz = (window - profile[i] / 2) * (ASigma[j] * ASigma[i])
      + window * AMean[j] * AMean[i];
    double lastzz = lastz;
    for (int k = 1; k < stepSize && i + k < timeSeriesLength - window
        + 1 && j + k < timeSeriesLength - window + 1; k++)
    {
      lastz = lastz - timeSeries_in[i + k - 1] * timeSeries_in[j + k - 1]
        + timeSeries_in[i + k + window - 1] * timeSeries_in[j
        + k + window - 1];
      distance = 2 * (window - (lastz - window * AMean[j + k] * AMean[i
            + k]) / (ASigma[j + k] * ASigma[i + k]));
      if (distance < profile[i + k])
      {
        profile[i + k] = distance;
        profileIndex[i + k] = j + k;
      }
      if (distance < profile[j + k])
      {
        profile[j + k] = distance;
        profileIndex[j + k] = i + k;
      }
    }
    lastz = lastzz;
    for (int k = 1; k < stepSize && i - k >= 0 && j - k >= 0; k++)
    {
      lastz = lastz - timeSeries_in[i - k + window] * timeSeries_in[j
        - k + window] + timeSeries_in[i - k] * timeSeries_in[j - k];
      distance = 2 * (window - (lastz - window * AMean[j - k] * AMean[i
            - k]) / (ASigma[j - k] * ASigma[i - k]));
      if (distance < profile[i - k])
      {
        profile[i - k] = distance;
        profileIndex[i - k] = j - k;
      }
      if (distance < profile[j - k])
      {
        profile[j - k] = distance;
        profileIndex[j - k] = i - k;
      }
    }
  }

  fftw_free(AFreq);
  fftw_free(queryTime);
  fftw_free(queryFreq);
  fftw_free(AQueryTime);
  fftw_free(AQueryFreq);

  /********************************** SCRIMP *********************************/

  //Random shuffle the computation order of the diagonals of the distance
  //matrix
  idx.clear();
  for (int i = exclusionZone+1; i < ProfileLength; i++)
    idx.push_back(i);
  shuffle(idx.begin(), idx.end(), g);

  vector<double> dotproduct(timeSeriesLength);

  //iteratively evaluate the diagonals of the distance matrix
  for (int ri = 0; ri < (int) idx.size(); ri++)
  {
    //select a random diagonal
    int diag = idx[ri];

    //calculate the dot product of every two time series values that ar diag
    //away
    for (int j=diag; j < timeSeriesLength; j++)
      dotproduct[j]=timeSeries_in[j]*timeSeries_in[j-diag];

    //evaluate the fist distance value in the current diagonal
    double distance;
    double lastz=0; //the dot product of a subsequence
    for (int k = 0; k < window; k++)
      lastz += dotproduct[k+diag];

    //j is the column index, i is the row index of the current distance value
    //in the distance matrix
      int j=diag, i=j-diag;

    //evaluate the distance based on the dot product
      distance = 2 * (window - (lastz - window * AMean[j] * AMean[i])
          / (ASigma[j] * ASigma[i]));

    //update matrix profile and matrix profile index if the current distance
    //value is smaller
    if (distance < profile[j])
    {
      profile[j] = distance;
      profileIndex [j] = i;
    }
    if (distance < profile[i])
    {
      profile[i] = distance;
      profileIndex [i] = j;
    }

    //evaluate the second to the last distance values along the diagonal and
    //update the matrix profile/matrix profile index.
    for (j=diag+1; j<ProfileLength; j++)
    {
      i=j-diag;
      lastz = lastz + dotproduct[j+window-1] - dotproduct [j-1];
      distance = 2 * (window - (lastz - window * AMean[j] * AMean[i])
          / (ASigma[j] * ASigma[i]));
      if (distance < profile[j])
      {
        profile[j] = distance;
        profileIndex [j] = i;
      }
      if (distance < profile[i])
      {
        profile[i] = distance;
        profileIndex [i] = j;
      }
    }
  }

  double ssf = numeric_limits<double>::infinity();

  //Write final Matrix Profile and Matrix Profile Index to file.
  for (int i = 0; i < ProfileLength; i++)
  {
    profile[i] = profile[i];
    if (profile[i] < ssf) {

      pos0_out = i;
      pos1_out = profileIndex[i];
      ssf = sqrt(profile[i]);
    }
  }

  fftw_cleanup();
}
