///\file unittests.hpp
///
///\brief File contains class definitions.
///
///This is the header file of the unit tests. This file contains class
///definitions with public functions to test private or protected functions of
///a class.

#ifndef UNITTESTS_HPP
#define UNITTESTS_HPP

#include <stest.hpp>
#include <tsgenerator.hpp>
#include <freepositions.hpp>
#include <basets.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <filesystem>

using namespace std;
using namespace std::filesystem;


// --------------------------------------------------------------------------
///\brief This variable contains the window size of the top motif pair and set.
///
///This variable stores the window size of the top motif pair and top motif set
///in the test time series.
// --------------------------------------------------------------------------
int windowSize = 20;

// --------------------------------------------------------------------------
///\brief This variable contains the top motif pair similarity.
///
///This variable stores the similarity of the subsequences of the top motif
///pair.
// --------------------------------------------------------------------------
double topMotifPairSimilarity = 0.0;

// --------------------------------------------------------------------------
///\brief This variable contains the positions of the top motif pair.
///
///This variable stores the positions of the top motif pair subsequences.
// --------------------------------------------------------------------------
vector<int> topMotifPairPos = { 150, 250 };


// --------------------------------------------------------------------------
///\brief This variable contains the range of the top motif set.
///
///This variable stores the range of the top motif set in the test time series.
// --------------------------------------------------------------------------
double topMotifSetRange = 0.240683084;

// --------------------------------------------------------------------------
///\brief This variable contains the top motif set size.
///
///This variable stores the  top motif set size, i. e. the number of
///subsequences in the top motif set.
// --------------------------------------------------------------------------
double topMotifSetSize = 3;

// --------------------------------------------------------------------------
///\brief This variable contains the positions of the top motif set.
///
///This variable stores the positions of the top motif set in the test time
///series.
// --------------------------------------------------------------------------
vector<int> topMotifSetPos = { 35, 150, 250 };

// --------------------------------------------------------------------------
///\brief This variable contains the means of the top motif set.
///
///This variable stores the means of the top motif set in the test time series.
// --------------------------------------------------------------------------
vector<double> topMotifSetMeans = { 31.7559967, 73.8733465, 73.8733465 };

// --------------------------------------------------------------------------
///\brief This variable contains the standard deviations of the top motif set.
///
///This variable stores the standard deviations of the top motif set in the
///test time series.
// --------------------------------------------------------------------------
vector<double> topMotifSetStdDevs = { 15.2633206, 16.6681750, 16.6681750 };

// --------------------------------------------------------------------------
///\brief This variable contains the test time series.
///
///This variable stores the test times series.
// --------------------------------------------------------------------------
vector<double> testTimeSeries({
    -10.768911,
    -18.433098,
    -11.754007,
    -5.639954,
    -14.184551,
    -12.368739,
    -22.100594,
    -21.991726,
    -24.945935,
    -24.239062,
    -15.301486,
    -27.441897,
    -23.450026,
    -31.510731,
    -39.267835,
    -31.694873,
    -28.583959,
    -29.351509,
    -25.496541,
    -34.101213,
    -17.759350,
    -19.176013,
    -18.355335,
    -4.800288,
    -10.434776,
    -16.000721,
    -24.525849,
    -22.024235,
    -27.125363,
    -13.327274,
    -17.874286,
    -16.360775,
    -21.200392,
    -16.082495,
    -24.664626,
    -16.603681,
    39.386230,
    36.565618,
    39.793586,
    43.013628,
    38.668908,
    33.500841,
    39.933202,
    35.754748,
    35.941286,
    37.664790,
    36.074232,
    38.735548,
    39.152116,
    34.042679,
    29.797315,
    32.463496,
    31.694503,
    38.783557,
    -9.242668,
    -20.720239,
    -16.059195,
    -16.283937,
    -27.705698,
    -27.528500,
    -17.071427,
    -29.126009,
    -32.409610,
    -30.973661,
    -37.291269,
    -36.608783,
    -26.005272,
    -37.764451,
    -32.923841,
    -27.450907,
    -19.964810,
    -30.979881,
    -28.806412,
    -31.286976,
    -29.794595,
    -32.534666,
    -41.688122,
    -35.074693,
    -35.577225,
    -40.667025,
    -42.325801,
    -40.884100,
    -27.198010,
    -27.086483,
    -27.380737,
    -32.583106,
    -30.999783,
    -32.115791,
    -34.640235,
    -29.296114,
    -16.688001,
    -27.592166,
    -26.411253,
    -23.397623,
    -15.070153,
    -32.235968,
    -15.874426,
    -8.023144,
    -12.147740,
    -11.254056,
    -10.229766,
    -10.044629,
    -14.556040,
    -7.515845,
    -4.536315,
    5.084161,
    6.216170,
    12.597693,
    10.858761,
    0.907416,
    5.198124,
    -6.700793,
    -6.081909,
    1.099458,
    -0.443738,
    5.349589,
    8.763887,
    -7.506228,
    -4.146017,
    -10.714550,
    4.991264,
    -0.873695,
    4.840859,
    7.151300,
    -8.191696,
    3.565129,
    4.786711,
    8.768263,
    4.033746,
    10.034495,
    11.365858,
    9.754581,
    18.323738,
    4.900006,
    15.448548,
    8.925523,
    15.782048,
    17.449776,
    16.828506,
    16.523008,
    17.702186,
    12.703089,
    16.764403,
    18.802697,
    26.338569,
    24.982530,
    22.969294,
    22.401425,
    27.509449,
    26.730151,
    19.884357,
    83.765535,
    81.204662,
    81.243895,
    85.012883,
    84.519502,
    71.486461,
    84.425170,
    77.610754,
    79.841677,
    81.416664,
    76.801356,
    82.202635,
    81.351707,
    77.429990,
    69.881711,
    74.637546,
    72.899982,
    80.076846,
    31.773598,
    27.382474,
    23.286007,
    33.989471,
    41.398836,
    26.797881,
    36.269419,
    23.224410,
    44.731778,
    32.596313,
    21.490972,
    22.313224,
    23.485559,
    16.545780,
    19.378278,
    17.164862,
    18.700601,
    21.415953,
    17.736355,
    22.207187,
    17.129972,
    15.667731,
    21.391168,
    16.843302,
    17.727531,
    8.930618,
    8.752489,
    11.471171,
    1.314668,
    7.748010,
    7.723902,
    12.359278,
    16.541422,
    7.933293,
    3.118250,
    -4.884358,
    -19.442805,
    -11.814778,
    -8.679566,
    -14.217084,
    -11.802162,
    -10.339839,
    -8.180208,
    -14.966384,
    -20.195693,
    -19.937768,
    -17.339444,
    -15.369381,
    -23.204906,
    -40.890050,
    -15.865825,
    -19.210405,
    -26.727909,
    -30.134980,
    -28.058161,
    -23.600434,
    -23.201012,
    -16.281554,
    -19.910757,
    -15.558630,
    -14.118548,
    -5.580549,
    -7.821060,
    -22.087571,
    -1.921322,
    -23.451956,
    -6.153180,
    -13.010591,
    -16.130146,
    -13.849764,
    -20.773420,
    -15.399393,
    -21.287319,
    -18.676833,
    -26.591504,
    -33.562469,
    -29.948043,
    -29.049934,
    -35.790035,
    -19.763629,
    -31.446928,
    19.884357,
    83.765535,
    81.204662,
    81.243895,
    85.012883,
    84.519502,
    71.486461,
    84.425170,
    77.610754,
    79.841677,
    81.416664,
    76.801356,
    82.202635,
    81.351707,
    77.429990,
    69.881711,
    74.637546,
    72.899982,
    80.076846,
    31.773598,
    -23.640862,
    -26.614654,
    -24.485532,
    -20.474911,
    -11.680886,
    -26.061593,
    -19.460513,
    -14.645193,
    16.773669,
    3.633842,
    11.245365,
    14.480009,
    17.406325,
    13.715258,
    13.006272,
    14.420967,
    6.342660,
    4.877567,
    9.314946,
    5.125607,
    11.778980,
    8.371439,
    6.399086,
    16.394781,
    4.545560,
    6.393353,
    8.814702,
    -0.494955,
    4.714512,
    0.685090 });


// --------------------------------------------------------------------------
///\brief This variable contains the test time series rolling means.
///
///This variable stores the test times series rolling means.
// --------------------------------------------------------------------------
vector<double> testMeans({
    -22.631332,
    -22.980854,
    -23.018000,
    -23.348066,
    -23.306083,
    -23.118594,
    -23.300194,
    -23.421456,
    -23.423082,
    -23.532053,
    -22.986464,
    -23.115104,
    -22.561048,
    -22.448566,
    -21.677154,
    -20.946994,
    -20.192434,
    -16.793925,
    -13.498068,
    -10.233562,
    -6.377820,
    -3.556407,
    -0.922564,
    1.991863,
    4.019614,
    6.338417,
    9.021693,
    12.051697,
    15.089686,
    18.403560,
    20.772058,
    23.155638,
    25.596851,
    28.241596,
    30.984899,
    31.755997,
    31.550169,
    28.777898,
    26.135420,
    22.760456,
    19.233349,
    16.446332,
    13.314990,
    9.697849,
    6.361429,
    2.699801,
    -1.013878,
    -4.117853,
    -7.942853,
    -11.546651,
    -14.621330,
    -17.109436,
    -20.281605,
    -23.306651,
    -26.810177,
    -27.837774,
    -28.428495,
    -29.709941,
    -30.649479,
    -31.043056,
    -31.699982,
    -32.962700,
    -33.550605,
    -33.290025,
    -33.095666,
    -32.600140,
    -32.398856,
    -32.648581,
    -32.366148,
    -32.451968,
    -32.544228,
    -32.380388,
    -32.211002,
    -32.091244,
    -31.696776,
    -30.960554,
    -30.945619,
    -29.654935,
    -28.302357,
    -27.130883,
    -25.660235,
    -24.055433,
    -22.513459,
    -21.881361,
    -20.902829,
    -19.760608,
    -17.877244,
    -16.016447,
    -13.780773,
    -11.505823,
    -9.995646,
    -8.901340,
    -7.856771,
    -6.840304,
    -5.615450,
    -4.884129,
    -3.004851,
    -1.772936,
    -1.747090,
    -1.347004,
    -1.320029,
    -0.558977,
    -0.100430,
    0.869415,
    1.602772,
    1.420003,
    1.344051,
    1.272578,
    1.081107,
    0.739856,
    1.196210,
    1.504597,
    2.327365,
    3.547648,
    3.737675,
    4.532289,
    4.711086,
    5.061994,
    6.309794,
    7.358520,
    8.720398,
    9.355944,
    10.034784,
    10.630961,
    11.213531,
    12.940044,
    14.010914,
    14.920043,
    15.601701,
    16.775486,
    17.610269,
    18.036194,
    21.736742,
    24.880788,
    28.697983,
    32.176199,
    35.955898,
    38.741119,
    42.089889,
    45.129001,
    48.294934,
    51.480658,
    54.685572,
    57.957483,
    61.084934,
    63.639505,
    65.884464,
    68.467877,
    70.992804,
    73.621174,
    73.873347,
    74.248252,
    71.224276,
    68.863516,
    66.871263,
    63.960513,
    61.548009,
    59.134907,
    57.150237,
    54.899515,
    51.981980,
    49.026808,
    46.361018,
    43.078175,
    39.979504,
    36.966247,
    34.407192,
    31.746112,
    28.987931,
    26.094448,
    25.362267,
    24.776529,
    24.681788,
    23.824479,
    22.640914,
    21.747551,
    20.371704,
    19.784042,
    17.613187,
    16.370772,
    15.682418,
    15.184721,
    14.837514,
    14.406890,
    13.593888,
    12.491427,
    10.584257,
    8.922720,
    7.601924,
    5.780711,
    4.334104,
    3.033726,
    1.555157,
    -0.035328,
    -1.931489,
    -3.374908,
    -4.679505,
    -6.021532,
    -7.247511,
    -9.679414,
    -10.858900,
    -12.437385,
    -14.600851,
    -16.504265,
    -18.063085,
    -18.998889,
    -19.186799,
    -19.410138,
    -19.971698,
    -20.038775,
    -20.154594,
    -19.916630,
    -19.898672,
    -20.254732,
    -19.341013,
    -19.516723,
    -18.957409,
    -18.839470,
    -18.485732,
    -17.133718,
    -17.379097,
    -17.188547,
    -16.916517,
    -16.343610,
    -16.270277,
    -16.768379,
    -17.105730,
    -17.744149,
    -18.538113,
    -18.748363,
    -19.614782,
    -18.341537,
    -13.762207,
    -8.597596,
    -4.439335,
    0.983907,
    5.517541,
    9.742394,
    14.770160,
    19.343186,
    24.373940,
    29.214743,
    34.119177,
    39.163150,
    44.560311,
    50.109934,
    55.101422,
    60.285796,
    65.720297,
    70.712320,
    73.873347,
    71.697086,
    66.178076,
    60.893566,
    55.807626,
    50.972938,
    45.443883,
    40.896534,
    35.943016,
    32.901162,
    29.090770,
    25.582205,
    22.466138,
    19.226322,
    15.844500,
    12.623314,
    9.850277,
    6.435532,
    3.034412,
    -0.503683,
    -1.836083,
    -0.065091,
    1.684214,
    3.228445,
    5.071929,
    5.883252,
    7.505999,
    8.919760,
    9.627272,
    9.024314,
    8.876876 });


// --------------------------------------------------------------------------
///\brief This variable contains the test time series running standard
///deviations.
///
///This variable stores the test times series running standard deviations.
// --------------------------------------------------------------------------
vector<double> testStds({
    8.620709,
    8.267129,
    8.248259,
    7.916311,
    8.011764,
    8.262966,
    8.062188,
    8.061472,
    8.061187,
    8.095693,
    8.391928,
    8.292267,
    8.354621,
    8.357039,
    8.195439,
    7.183783,
    6.797408,
    14.443506,
    18.227166,
    21.362922,
    23.553930,
    25.333982,
    26.293286,
    27.406343,
    28.314057,
    28.927706,
    29.218614,
    28.720491,
    28.163454,
    26.870856,
    26.044576,
    24.536386,
    22.854522,
    20.191446,
    17.535512,
    15.263321,
    15.927380,
    18.874807,
    21.160616,
    23.916465,
    25.797213,
    26.547006,
    28.004238,
    28.987310,
    29.629258,
    30.265738,
    30.304360,
    29.515617,
    28.658768,
    26.993394,
    25.058284,
    22.902047,
    20.029648,
    16.142915,
    7.664438,
    6.534686,
    6.397128,
    6.357833,
    5.653801,
    5.708888,
    6.014423,
    5.433645,
    5.619634,
    5.784899,
    5.923108,
    5.965780,
    5.894622,
    5.721728,
    5.600357,
    5.621357,
    5.553249,
    5.955724,
    6.040715,
    6.130105,
    6.416322,
    7.366715,
    7.363810,
    7.625406,
    8.845642,
    9.342132,
    9.410363,
    9.165013,
    8.790603,
    8.885040,
    9.324706,
    9.845834,
    10.771987,
    11.532030,
    12.488734,
    12.625076,
    12.206212,
    12.533845,
    11.780531,
    10.985981,
    10.422872,
    10.245447,
    8.322728,
    8.148254,
    8.129185,
    7.797606,
    7.764148,
    7.597713,
    7.281555,
    6.546223,
    6.385350,
    6.606924,
    6.573010,
    6.527256,
    6.242176,
    5.874018,
    6.214011,
    6.548981,
    6.499915,
    7.072355,
    7.055059,
    7.424653,
    7.484993,
    7.823621,
    7.708883,
    7.641532,
    6.663776,
    6.880428,
    6.496719,
    6.539685,
    6.720235,
    5.898381,
    6.041586,
    5.952527,
    5.989512,
    5.907227,
    6.072988,
    5.916816,
    15.293878,
    20.006422,
    22.903370,
    25.734362,
    27.530577,
    28.159315,
    29.383981,
    29.755031,
    29.911069,
    29.876075,
    28.968609,
    28.185698,
    27.117563,
    26.111855,
    24.576804,
    22.563034,
    19.939697,
    17.328274,
    16.668175,
    15.492483,
    18.873293,
    20.370791,
    21.001233,
    22.280422,
    22.534491,
    23.884649,
    23.343740,
    23.432404,
    23.775334,
    23.605618,
    23.327350,
    22.662845,
    21.420641,
    20.141068,
    19.016418,
    16.794563,
    14.127328,
    7.937520,
    8.054414,
    8.308167,
    8.335383,
    8.214869,
    7.245618,
    7.761152,
    7.499458,
    7.710422,
    6.377649,
    5.724700,
    5.892873,
    5.729925,
    5.418332,
    5.604497,
    5.990438,
    7.148749,
    9.824937,
    10.629592,
    11.083851,
    11.518364,
    11.815143,
    11.926873,
    11.379893,
    11.354912,
    11.396695,
    11.752159,
    11.781680,
    11.387643,
    11.842627,
    13.404529,
    12.847569,
    11.794134,
    10.131308,
    9.257240,
    8.407678,
    7.915964,
    7.968702,
    7.820168,
    7.422597,
    7.376207,
    7.263290,
    7.648687,
    7.676591,
    7.604372,
    8.590524,
    8.636748,
    9.108939,
    9.169700,
    9.130851,
    7.584289,
    7.618610,
    7.618083,
    7.365542,
    6.733612,
    6.612533,
    7.466051,
    7.889509,
    8.302771,
    9.184423,
    9.161913,
    9.496361,
    12.518819,
    25.524680,
    32.746080,
    38.162294,
    42.531872,
    46.203562,
    48.138863,
    50.373387,
    51.701620,
    52.443163,
    53.013674,
    52.650713,
    52.181279,
    50.661293,
    47.797810,
    44.258308,
    39.961614,
    33.374008,
    27.089381,
    16.668175,
    24.552035,
    32.377778,
    37.684202,
    41.286401,
    43.200467,
    45.564868,
    47.246133,
    47.614775,
    46.791849,
    45.908781,
    44.433424,
    42.890700,
    40.644475,
    38.066745,
    35.347775,
    32.833075,
    29.276202,
    24.995471,
    17.817151,
    16.284007,
    15.733040,
    14.587093,
    13.314179,
    12.427553,
    11.822303,
    9.280242,
    6.917382,
    4.900622,
    4.722896,
    4.930343 });


// --------------------------------------------------------------------------
///\brief This variable contains the larger motif set test time series.
///
///This variable stores the test times series to test the function that checks
///if there is a larger motif set.
// --------------------------------------------------------------------------
vector<double> testTimeSeriesCheckLargerMotifSet( {
    -10.768911,
    -18.433098,
    -11.754007,
    -5.639954,
    -14.184551,
    -12.368739,
    -22.100594,
    -21.991726,
    -24.945935,
    -24.239062,
    -15.301486,
    -27.441897,
    -23.450026,
    -31.510731,
    -39.267835,
    -31.694873,
    -28.583959,
    -29.351509,
    -25.496541,
    -34.101213,
    -17.759350,
    -19.176013,
    -18.355335,
    -4.800288,
    -10.434776,
    -16.000721,
    -24.525849,
    -22.024235,
    -27.125363,
    -13.327274,
    -17.874286,
    -16.360775,
    -21.200392,
    -16.082495,
    -24.664626,
    -16.603681,
    39.386230,
    36.565618,
    39.793586,
    43.013628,
    38.668908,
    33.500841,
    39.933202,
    35.754748,
    35.941286,
    37.664790,
    36.074232,
    38.735548,
    39.152116,
    34.042679,
    29.797315,
    32.463496,
    31.694503,
    38.783557,
    -9.242668,
    -20.720239,
    -16.059195,
    -16.283937,
    -27.705698,
    -27.528500,
    -17.071427,
    -29.126009,
    -32.409610,
    -30.973661,
    -37.291269,
    -36.608783,
    -26.005272,
    -37.764451,
    -32.923841,
    -27.450907,
    -19.964810,
    -30.979881,
    -28.806412,
    -31.286976,
    -29.794595,
    -32.534666,
    -41.688122,
    -35.074693,
    -35.577225,
    -40.667025,
    -42.325801,
    -40.884100,
    -27.198010,
    -27.086483,
    -27.380737,
    -32.583106,
    -30.999783,
    -32.115791,
    -34.640235,
    -29.296114,
    -16.688001,
    -27.592166,
    -26.411253,
    -23.397623,
    -15.070153,
    -32.235968,
    -15.874426,
    -8.023144,
    -12.147740,
    -11.254056,
    -10.229766,
    -10.044629,
    -14.556040,
    -7.515845,
    -4.536315,
    5.084161,
    6.216170,
    12.597693,
    10.858761,
    0.907416,
    5.198124,
    -6.700793,
    -6.081909,
    1.099458,
    -0.443738,
    5.349589,
    8.763887,
    -7.506228,
    -4.146017,
    -10.714550,
    4.991264,
    -0.873695,
    4.840859,
    7.151300,
    -8.191696,
    3.565129,
    4.786711,
    8.768263,
    4.033746,
    19.884357,
    83.765535,
    81.204662,
    81.243895,
    85.012883,
    84.519502,
    71.486461,
    84.425170,
    77.610754,
    79.841677,
    81.416664,
    76.801356,
    82.202635,
    81.351707,
    77.429990,
    69.881711,
    74.637546,
    72.899982,
    80.076846,
    31.773598,
    26.730151,
    19.884357,
    83.765535,
    81.204662,
    81.243895,
    85.012883,
    84.519502,
    71.486461,
    84.425170,
    77.610754,
    79.841677,
    81.416664,
    76.801356,
    82.202635,
    81.351707,
    77.429990,
    69.881711,
    74.637546,
    72.899982,
    80.076846,
    31.773598,
    27.382474,
    23.286007,
    33.989471,
    41.398836,
    26.797881,
    36.269419,
    23.224410,
    44.731778,
    32.596313,
    21.490972,
    22.313224,
    23.485559,
    16.545780,
    19.378278,
    17.164862,
    18.700601,
    21.415953,
    17.736355,
    22.207187,
    17.129972,
    15.667731,
    21.391168,
    16.843302,
    17.727531,
    8.930618,
    8.752489,
    11.471171,
    1.314668,
    7.748010,
    7.723902,
    12.359278,
    16.541422,
    7.933293,
    3.118250,
    -4.884358,
    -19.442805,
    -11.814778,
    -8.679566,
    -14.217084,
    -11.802162,
    -10.339839,
    -8.180208,
    -14.966384,
    -20.195693,
    -19.937768,
    -17.339444,
    -15.369381,
    -23.204906,
    -40.890050,
    -15.865825,
    -19.210405,
    -26.727909,
    -30.134980,
    -28.058161,
    -23.600434,
    -23.201012,
    -16.281554,
    -19.910757,
    -15.558630,
    19.884357,
    83.765535,
    81.204662,
    81.243895,
    85.012883,
    84.519502,
    71.486461,
    84.425170,
    77.610754,
    79.841677,
    81.416664,
    76.801356,
    82.202635,
    81.351707,
    77.429990,
    69.881711,
    74.637546,
    72.899982,
    80.076846,
    31.773598,
    -31.446928,
    19.884357,
    83.765535,
    81.204662,
    81.243895,
    85.012883,
    84.519502,
    71.486461,
    84.425170,
    77.610754,
    79.841677,
    81.416664,
    76.801356,
    82.202635,
    81.351707,
    77.429990,
    69.881711,
    74.637546,
    72.899982,
    80.076846,
    31.773598,
    -23.640862,
    -26.614654,
    -24.485532,
    -20.474911,
    -11.680886,
    -26.061593,
    -19.460513,
    -14.645193,
    16.773669,
    3.633842,
    11.245365,
    14.480009,
    17.406325,
    13.715258,
    13.006272,
    14.420967,
    6.342660,
    4.877567,
    9.314946,
    5.125607,
    11.778980,
    8.371439,
    6.399086,
    16.394781,
    4.545560,
    6.393353,
    8.814702,
    -0.494955,
    4.714512,
    0.685090 });

// --------------------------------------------------------------------------
///\brief This variable contains the test top motif set subsequences.
///
///This variable stores the test top motif set subsequences.
// --------------------------------------------------------------------------
vector<vector<double>> testMotifSetSubsequences = {
  {
    -16.603681,
    39.386230,
    36.565618,
    39.793586,
    43.013628,
    38.668908,
    33.500841,
    39.933202,
    35.754748,
    35.941286,
    37.664790,
    36.074232,
    38.735548,
    39.152116,
    34.042679,
    29.797315,
    32.463496,
    31.694503,
    38.783557,
    -9.242668
  },
  {
    19.884357,
    83.765535,
    81.204662,
    81.243895,
    85.012883,
    84.519502,
    71.486461,
    84.425170,
    77.610754,
    79.841677,
    81.416664,
    76.801356,
    82.202635,
    81.351707,
    77.429990,
    69.881711,
    74.637546,
    72.899982,
    80.076846,
    31.773598
  },
  {
    19.884357,
    83.765535,
    81.204662,
    81.243895,
    85.012883,
    84.519502,
    71.486461,
    84.425170,
    77.610754,
    79.841677,
    81.416664,
    76.801356,
    82.202635,
    81.351707,
    77.429990,
    69.881711,
    74.637546,
    72.899982,
    80.076846,
    31.773598
  }
};


// --------------------------------------------------------------------------
///\brief This variable contains the top motif pair similarity.
///
///This variable stores the similarity of the subsequences of the top motif
///pair.
// --------------------------------------------------------------------------
double testSequencesSimilarty = 0.4813661;

// --------------------------------------------------------------------------
///\brief This variable contains the first test sequence.
///
///This variable stores the first test sequence.
// --------------------------------------------------------------------------
vector<double> testSequenceOne = {
  -3.1683588980060402,
  0.49990650723292546,
  0.31510976016894321,
  0.52659506407697054,
  0.73756108386608021,
  0.45291004965654513,
  0.11431616641109528,
  0.53574222172860275,
  0.26198436103862266,
  0.27420568591557165,
  0.38712370964650772,
  0.2829158482261514,
  0.45727607207449689,
  0.48456816731161162,
  0.14981551897559328,
  -0.12832605359891772,
  0.046352908230570186,
  -0.0040288546332950858,
  0.46042145578193161,
  -2.6860907741039628 };

// --------------------------------------------------------------------------
///\brief This variable contains the second test sequence.
///
///This variable stores the second test sequence.
// --------------------------------------------------------------------------
vector<double> testSequenceTwo = {
  -3.2390462337465022,
  0.5934775962570914,
  0.43983912080329196,
  0.44219288778333282,
  0.66831170364170056,
  0.63871152250898888,
  -0.14320017313810549,
  0.63305211469512646,
  0.22422415432849474,
  0.3580674212827003,
  0.45255808971473588,
  0.17566466904540784,
  0.49971208219733798,
  0.44866102240602851,
  0.21337929580683826,
  -0.23947645996860673,
  0.045847817693667449,
  -0.058396588008374987,
  0.37217628440632827,
  -2.5257563277094572 };



// --------------------------------------------------------------------------
///\brief This class represents a test version of the time series generator.
///
///The protected functions as well as the protected attributes of the test time
///series generator are made accessable.
// --------------------------------------------------------------------------
class TestTSGenerator : public tsg::TSGenerator {

public:
  // --------------------------------------------------------------------------
  ///\brief The constructor initializes the TSGenerator.
  ///
  ///\param [in] length_in Hands over the time series length.
  ///\param [in] window_in Hands over the window size.
  ///\param [in] delta_in Hands over the maximum difference between two
  ///consecutive values in the time series.
  ///\param [in] noise_in Hands over the noise option, i.e., +-noise_in 2 will
  ///be added to each value of the time series.
  ///\param [in] type_in Hands over the motif type, i.e., the motif shape.
  ///\param [in] size_in Hands over the number of subsequences non-self matched
  ///by the motif.
  ///\param [in] height_in Hands over the maximum difference between two values
  ///of the motif.
  ///
  ///The constructor checks whether a true random engine is available and
  ///stores the result in the trueRandomEngineAvailable variable.
  // --------------------------------------------------------------------------
  TestTSGenerator(int length_in, int window_in, double delta_in, double
      noise_in, int type_in, int size_in, double height_in)
    : TSGenerator(length_in, window_in, delta_in, noise_in, type_in, size_in,
        height_in) { }

  // --------------------------------------------------------------------------
  ///\brief Sets the length of the time series.
  ///
  ///\param [in] length_in The legnth of the time series.
  ///
  ///This function sets the content of the variable that stores the length
  ///of the time series.
  // --------------------------------------------------------------------------
  void setLength(int length_in) {

    length = length_in;
  }

  // --------------------------------------------------------------------------
  ///\brief Sets the window size for the motif detection.
  ///
  ///\param [in] window_in The window size.
  ///
  ///This function sets the content of the variable that stores the window
  ///size for the motif detection.
  // --------------------------------------------------------------------------
  void setWindow(int window_in) {

    window = window_in;
  }

  // --------------------------------------------------------------------------
  ///\brief Returns the running sum of the time series.
  ///
  ///\return The running sums of the time series.
  ///
  ///This function returns the content of the variable that stores the sums of
  ///the time series.
  // --------------------------------------------------------------------------
  vector<double> &getSums() {

    return sums;
  }

  // --------------------------------------------------------------------------
  ///\brief Returns the running sum of square of the time series.
  ///
  ///\return The sum of squares of the time series.
  ///
  ///This function returns the content of the variable that stores the
  ///sum of squares of the time series.
  // --------------------------------------------------------------------------
  vector<double> &getSumSquares() {

    return sumSquares;
  }

  // --------------------------------------------------------------------------
  ///\brief Runs the running sum and sum of square function.
  ///
  ///\param [in] &sequence_in Hands over the sequence.
  ///
  ///This function runs the running sum and sum of square function since the
  ///function is protected.
  // --------------------------------------------------------------------------
  void testCalcRunnings(const vector<double> &sequence_in) {

    calcRunnings(sequence_in);
  }

  // --------------------------------------------------------------------------
  ///\brief Runs the mean and variance function.
  ///
  ///\param [in] &sequence_in Hands over the sequence.
  ///\param [in] pos_in Hands over the position of changed subsequence.
  ///
  ///This function runs the update running sum and sum of square function since
  ///the function is protected.
  // --------------------------------------------------------------------------
  void testUpdateRunnings(const vector<double> &sequence_in, const int pos_in)
    {

    updateRunnings(sequence_in, pos_in);
  }

  // --------------------------------------------------------------------------
  ///\brief Runs the similarity function.
  ///
  ///\param [in] &timeSeriesOne_in Hands over the time series.
  ///\param [in] subsequenceOnePos_in Hands over the position of the first
  ///subsequence in the time series.
  ///\param [in] subsequenceTwoPos_in Hands over the position of the second
  ///subsequence in the time series.
  ///\param [in] bestSoFar_in Hands over the best similarity so far.
  ///
  ///\return The similarity of the two z-normalized time series.
  ///
  ///This function runs the similarity function since the similarity function
  ///is protected.
  // --------------------------------------------------------------------------
  double testSimilarity(const vector<double> &timeSeries_in, const int
      subsequenceOnePos_in, const int subsequenceTwoPos_in, const double
      bestSoFar_in) {

    return similarity(timeSeries_in, subsequenceOnePos_in,
        subsequenceTwoPos_in, bestSoFar_in);
  }

  // --------------------------------------------------------------------------
  ///\brief Runs the mean and standard deviation function.
  ///
  ///\param [in] &sequence_in Hands over the sequence.
  ///\param [out] &mean_out The mean of the sequence.
  ///\param [out] &stdDev_out The standard deviation of the sequence.
  ///
  ///This function runs the mean and standard deviation function since the mean
  ///and standard deviation function is protected.
  // --------------------------------------------------------------------------
  void testMeanStdDev(const vector<double> &sequence_in, double &mean_out,
      double &stdDev_out) {

    meanStdDev(sequence_in, mean_out, stdDev_out);
  }

  // --------------------------------------------------------------------------
  ///\brief Runs the similarity function.
  ///
  ///\param [in] &timeSeriesOne_in Hands over the time series.
  ///\param [in] &subsequenceOne_in Hands over the first subsequence in the
  ///time series.
  ///\param [in] meanOne_in Hands over the mean of timeSeriesOne_in.
  ///\param [in] stdDevOne_in Hands over the standard deviation of
  ///timeSeriesOne_in.
  ///\param [in] subsequenceTwoPos_in Hands over the position of the second
  ///subsequence in the time series.
  ///\param [in] bestSoFar_in Hands over the best similarity so far.
  ///
  ///\return The similarity of the two z-normalized time series.
  ///
  ///This function runs the similarity function since the similarity function
  ///is protected.
  // --------------------------------------------------------------------------
  double testSimilarity(const vector<double> &timeSeries_in, const
      vector<double> &subsequenceOne_in, const double meanOne_in, const double
      stdDevOne_in, const int subsequenceTwoPos_in, const double bestSoFar_in)
    {

    return similarity(timeSeries_in, subsequenceOne_in, meanOne_in,
        stdDevOne_in, subsequenceTwoPos_in, bestSoFar_in);
  }

  // --------------------------------------------------------------------------
  ///\brief Runs the calculate motif set subsequence function.
  ///
  ///\param [out] subsequence_out Hands over the calculated raw subsequence.
  ///\param [in] type_in Hands over the motif set type.
  ///\param [in] height_in Hands over the height of the subsequence.
  ///
  ///This function runs the calculate motif set subsequence function since the
  ///calculate raw motif set subsequence function is protected.
  // --------------------------------------------------------------------------
  void testCalculateSubsequence(vector<double> &subsequence_out, int type_in,
      double height_in) {

    calculateSubsequence(subsequence_out, type_in, height_in);
  }

  // --------------------------------------------------------------------------
  ///\brief Runs the search for unintentional matches in the time series
  ///function.
  ///
  ///\param [in] &timeSereis_in Hands over the time series.
  ///\param [in] &subsequencePositions_in Hands over the position of the new
  ///subsequence.
  ///\param [in] similarity_in Hands over the similarity to break.
  ///
  ///\return true if there exists an unintentional subsequence match with the
  ///new subequence.
  ///
  ///This function runs the search for unintentional matches in the time series
  ///function since the function is protected.
  // --------------------------------------------------------------------------
  bool testSearchForUnintentionalMatches(const vector<double> &timeSeries_in,
      const vector<int> &motifSetPositions_in, double similarity_in) {

    return searchForUnintentionalMatches(timeSeries_in, motifSetPositions_in,
        similarity_in);
  }

  // --------------------------------------------------------------------------
  ///\brief Runs the larger motif set check function.
  ///
  ///\param [in] &timeSereis_in Hands over the time series.
  ///\param [in] &subsequencePositions_in Hands over the position of the new
  ///subsequence.
  ///\param [in] range_in Hands over the motif set range.
  ///
  ///\return true if there exists an unintentional subsequence match with the
  ///new subequence.
  ///
  ///This function runs the larger motif set in the set of motif set positions
  ///including all overlapping subsequences in the time series check function
  ///since the function is protected.
  // --------------------------------------------------------------------------
  bool testCheckIfThereIsALargerMotifSet(const vector<double> &timeSeries_in,
      const vector<int> &motifSetPositions_in,  double range_in) {

    return checkIfThereIsALargerMotifSet(timeSeries_in, motifSetPositions_in,
        range_in);
  }
};

#endif