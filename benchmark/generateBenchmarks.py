#!/usr/bin/python

import subprocess
import glob
import os
import shutil

lsb = "timeSeriesMotifBenchmark"
ldb = "timeSeriesMotifDatabaseBenchmark"

tsf = "time_series_"

def main():
  for path in glob.glob(tsf + "*"):
    if os.path.isdir(path):
      shutil.rmtree(path)
    else:
      os.remove(path)
  
  print("Generating motif benchmark.")
  print("***************************")
  
  for itrType in [ "box",
                   "triangle",
                   "semicircle",
                   "trapezoid",
                   "positiveflank",
                   "negativeflank",
                   "sine",
                   "cosine"
                 ]:
    for itrHeight in [ "150.0",
                       "200.0"
                     ]:
      for itrRandom in [ "0.03",
                         "0.04"
                       ]:
        for itrSize in [ "7",
                         "8"
                       ]:
          for itrLength in [ "100",
                             "150"
                           ]:
            print( "Generating time series with 10000 values, motif type "
                   + itrType
                   + " with "
                   + itrSize
                   + " non-self matching subsequences, subsequence length "
                   + itrLength
                   + ", subsequence height "
                   + itrHeight
                   + " and randomness factor "
                   + itrRandom
                   + "."
                 )
            while subprocess.call( [ "./TSGenerator",
                                     "-rd",
                                     itrRandom,
                                     "-l",
                                     "10000",
                                     "-w",
                                     itrLength,
                                     "-lm",
                                     itrType,
                                     itrSize,
                                     itrHeight
                                   ]
                                 ) != 0:
              print("retry")
  
  for itrRandom in [ "0.01",
                     "0.02",
                     "0.03",
                     "0.04",
                     "0.05",
                     "0.06",
                     "0.07",
                     "0.08"
                   ]:
    for itrType in [ "box",
                     "semicircle"
                   ]:
      for itrHeight in [ "150.0",
                         "200.0"
                       ]:
        for itrSize in [ "7",
                         "8"
                       ]:
          for itrLength in [ "100",
                             "150"
                           ]:
            print( "Generating time series with 10000 values, motif type "
                   + itrType
                   + " with "
                   + itrSize
                   + " non-self matching subsequences, subsequence length "
                   + itrLength
                   + ", subsequence height "
                   + itrHeight
                   + " and randomness factor "
                   + itrRandom
                   + "."
                 )
            while subprocess.call( [ "./TSGenerator",
                                     "-rd",
                                     itrRandom,
                                     "-l",
                                     "10000",
                                     "-w",
                                     itrLength,
                                     "-lm",
                                     itrType,
                                     itrSize,
                                     itrHeight,
                                   ]
                                 ) != 0:
              print("retry")
   
  for itrSize in [ "3",
                   "4",
                   "5",
                   "6",
                   "7",
                   "8",
                   "9",
                   "10"
                 ]:
    for itrType in [ "box",
                     "semicircle"
                   ]:
      for itrRandom in [ "0.02",
                         "0.03"
                       ]:
        for itrHeight in [ "150.0",
                           "200.0"
                         ]:
          for itrLength in [ "100",
                             "150"
                           ]:
            print( "Generating time series with 10000 values, motif type "
                   + itrType
                   + " with "
                   + itrSize
                   + " non-self matching subsequences, subsequence length "
                   + itrLength
                   + ", subsequence height "
                   + itrHeight
                   + " and randomness factor "
                   + itrRandom
                   + "."
                 )
            while subprocess.call( [ "./TSGenerator",
                                     "-rd",
                                     itrRandom,
                                     "-l",
                                     "10000",
                                     "-w",
                                     itrLength,
                                     "-lm",
                                     itrType,
                                     itrSize,
                                     itrHeight
                                   ]
                                 ) != 0:
              print("retry")
  
  if os.path.isdir(lsb):
    shutil.rmtree(lsb)
  os.mkdir(lsb)
  for path in glob.glob(tsf + "*"):
    shutil.move(path, lsb)
  
  print("Generating motif benchmark (database version).")
  print("**********************************************")
  
  for itrType in [ "box",
                   "triangle",
                   "semicircle",
                   "trapezoid",
                   "positiveflank",
                   "negativeflank",
                   "sine",
                   "cosine"
                 ]:
    for itrHeight in [ "150.0",
                       "200.0"
                     ]:
      for itrRandom in [ "0.03",
                         "0.04"
                       ]:
        for itrSize in [ "7",
                         "8"
                       ]:
          for itrLength in [ "100",
                             "150"
                           ]:
            print( "Generating time series with 10000 values, motif type "
                   + itrType
                   + " with "
                   + itrSize
                   + " non-self matching subsequences, subsequence length "
                   + itrLength
                   + ", subsequence height "
                   + itrHeight
                   + " and randomness factor "
                   + itrRandom
                   + "."
                 )
            while subprocess.call( [ "./TSGenerator",
                                     "-rd",
                                     itrRandom,
                                     "-l",
                                     "10000",
                                     "-w",
                                     itrLength,
                                     "-lm",
                                     itrType,
                                     itrSize,
                                     itrHeight
                                   ]
                                 ) != 0:
              print("retry")
  
  for itrRandom in [ "0.01",
                     "0.02",
                     "0.03",
                     "0.04",
                     "0.05",
                     "0.06",
                     "0.07",
                     "0.08"
                   ]:
    for itrType in [ "box",
                     "semicircle"
                   ]:
      for itrHeight in [ "150.0",
                         "200.0"
                       ]:
        for itrSize in [ "7",
                         "8"
                       ]:
          for itrLength in [ "100",
                             "150"
                           ]:
            print( "Generating time series with 10000 values, motif type "
                   + itrType
                   + " with "
                   + itrSize
                   + " non-self matching subsequences, subsequence length "
                   + itrLength
                   + ", subsequence height "
                   + itrHeight
                   + " and randomness factor "
                   + itrRandom
                   + "."
                 )
            while subprocess.call( [ "./TSGenerator",
                                     "-rd",
                                     itrRandom,
                                     "-l",
                                     "10000",
                                     "-w",
                                     itrLength,
                                     "-lm",
                                     itrType,
                                     itrSize,
                                     itrHeight,
                                     "-db"
                                   ]
                                 ) != 0:
              print("retry")
  
  for itrSize in [ "3",
                   "4",
                   "5",
                   "6",
                   "7",
                   "8",
                   "9",
                   "10"
                 ]:
    for itrType in [ "box",
                     "semicircle"
                   ]:
      for itrRandom in [ "0.02",
                         "0.03"
                       ]:
        for itrHeight in [ "150.0",
                           "200.0"
                         ]:
          for itrLength in [ "100",
                             "150"
                           ]:
            print( "Generating time series with 10000 values, motif type "
                   + itrType
                   + " with "
                   + itrSize
                   + " non-self matching subsequences, subsequence length "
                   + itrLength
                   + ", subsequence height "
                   + itrHeight
                   + " and randomness factor "
                   + itrRandom
                   + "."
                 )
            while subprocess.call( [ "./TSGenerator",
                                     "-rd",
                                     itrRandom,
                                     "-l",
                                     "10000",
                                     "-w",
                                     itrLength,
                                     "-lm",
                                     itrType,
                                     itrSize,
                                     itrHeight,
                                     "-db"
                                   ]
                                 ) != 0:
              print("retry")
  
  if os.path.isdir(ldb):
    shutil.rmtree(ldb)
  os.mkdir(ldb)
  for path in glob.glob(tsf + "*"):
    shutil.move(path, ldb)

main()
