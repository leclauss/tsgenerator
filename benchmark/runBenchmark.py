#!/usr/bin/python

import os
import sys
import subprocess
import glob
import shutil
from pathlib import Path
import numpy
import time

time_series_path = "timeSeriesMotifBenchmark/time_series_"

def getSubsequences( rs,
                     ofp,
                     ws
                   ):
  out = []
  
  sub = set()

  for line in rs:
    
    line = line.decode()

    ofp.write(line)

    if line == "\n" or line == "next\n" or line == "end\n":
      out += [sub]
      sub = set()
    else:
      pos = int(line)
      sub = sub.union(set(range(pos, pos + ws)))
  
  return out

def calcStats( injSet,
               calSet
             ):
  tp = len(injSet.intersection(calSet))
  fp = len(calSet.difference(injSet))
  fn = len(injSet.difference(calSet))

  return tp, fp, fn

def bestStats( injSet,
               calSets
             ):
  tp   = 0
  fp   = 0
  fn   = 0
  fscore = -1
  
  for calSet in calSets:
    ttp, tfp, tfn = calcStats(injSet, calSet)
    
    if ttp + tfp == 0:
      prec = 0
    else:
      prec = ttp / (ttp + tfp)
    
    if ttp + tfn == 0:
      rec = 0
    else:
      rec = ttp / (ttp + tfn)
    
    if prec + rec == 0:
      tfscore = 0
    else:
      tfscore = 2*prec*rec/(prec+rec)
    
    if fscore < tfscore:
      fscore = tfscore
      tp = ttp
      fp = tfp
      fn = tfn
 
  return tp, fp, fn

def save_output_and_create_stats( inj,
                                  rs,
                                  sfp,
                                  ofp,
                                  ws,
                                  tfp=-1,
                                  tb=-1
                                ):
  sub = getSubsequences(rs, ofp, ws)
  tp, tn, fn = bestStats(inj, sub)
  if tfp != -1 and tb != -1:
    tfp.write(str(time.time() - tb) + "\n")
  sfp.write(str(tp) + ", " + str(tn) + ", " + str(fn) + "\n")

  return tp, tn, fn

def main():
  path = sys.argv[0]

  my_file = Path(path)
  
  if path[0] == '~':
    os.chdir(os.path.dirname(os.path.join(str(Path.home()), path[1:])))
  elif my_file.is_file():
    os.chdir(os.path.dirname(os.path.join("./", path)))
  else:
    os.chdir(os.path.dirname(path))
  
  for path in glob.glob("*Stats.csv"):
    if os.path.isdir(path):
      os.remove(path)
  if os.path.isdir("stats.pdf"):
    os.remove("stats.pdf")
  for path in glob.glob("*Runtimes.csv"):
    if os.path.isdir(path):
      os.remove(path)
  
  #run motif discovery algorithms
  for i in range(0, 384):
    print(time_series_path + str(i))
    ts_abs_path = os.path.join( os.getcwd(),
                                time_series_path + str(i) + "/time_series_" + str(i) + ".csv"
                              )
    tsm_abs_path = os.path.join( os.getcwd(),
                                 time_series_path + str(i) + "/time_series_meta_" + str(i) + ".csv"
                               )
    mfp = open(tsm_abs_path, "r")
    mfp.readline()
    r = float(mfp.readline().split(", ")[-2])
    ws = int(mfp.readline().split(", ")[-2])
  
    #get the inj subsequences
    inj = set()
    for line in mfp:
      pos = int(line.split(", ")[-2])
      inj = inj.union(set(range(pos, pos + ws)))
  
    tb = time.time()
    #www.cs.ucr.edu/%7Emueen/zip/MK_code.zip
    #MK algorithm: mk_l file timeserieslength subsequencelength windowsize numberofreferencepoints
    #./mk timeSeriesTopPairMotifBenchmark/time_series_0.csv 10000 40 40 10
    rs = subprocess.Popen( [ "./algorithms/mk",
                             ts_abs_path,
                             "10000",
                             str(ws),
                             str(ws),
                             "10",
                             str(r)
                           ],
                           stdout = subprocess.PIPE
                         )
    sfp = open("mkStats.csv", "a+")
    tfp = open("mkRuntimes.csv", "a+")
    ofp = open(time_series_path + str(i) + "/mk.csv", "w")
    save_output_and_create_stats( inj,
                                  rs.stdout,
                                  sfp,
                                  ofp,
                                  ws,
                                  tfp,
                                  tb
                                )
  
    tb = time.time()
    ##https://github.com/javidlakha/matrix-profile
    #Matrix Profile I: mp -W ignore file windowsize
    #python -W ignore mp.py timeSeriesTopPairMotifBenchmark/time_series_0.csv 40
    rs = subprocess.Popen( [ "python",
                             "-W",
                             "ignore",
                             "algorithms/mp.py",
                             ts_abs_path,
                             str(ws),
                             str(r)
                           ],
                           stdout = subprocess.PIPE
                         )
    sfp = open("mpStats.csv", "a+")
    tfp = open("mpRuntimes.csv", "a+")
    ofp = open(time_series_path + str(i) + "/mp.csv", "w")
    save_output_and_create_stats( inj,
                                  rs.stdout,
                                  sfp,
                                  ofp,
                                  ws,
                                  tfp,
                                  tb
                                )
  
    tb = time.time()
    #http://fs.ismll.de/publicspace/LearnMotifs/
    #Learn Motifs: java -Xmx7g -jar LearnMotifs.jar dataSet=timeSeriesTopLatentMotifBenchmark/time_series_0.csv eta=0.3 maxIter=300 numRandomRestarts=1 alpha=2 tsLength=10000 w=200 K=3 pct=0.01
    #java -Xmx7g -jar LearnMotifs.jar dataSet=timeSeriesTopLatentMotifBenchmark/time_series_0.csv eta=0.3 maxIter=300 numRandomRestarts=1 alpha=2 tsLength=10000 w=200 K=3 pct=0.01
    rs = subprocess.Popen( [ "java",
                             "-Xmx7g",
                             "-jar",
                             "algorithms/LearnMotifs.jar",
                             "dataSet=" + ts_abs_path,
                             "eta=0.1",
                             "maxIter=1000",
                             "numRandomRestarts=200",
                             "alpha=2",
                             "tsLength=10000",
                             "w=" + str(ws),
                             "K=3",
                             "t=" + str(r*r)
                           ],
                           stdout = subprocess.PIPE
                         )
    sfp = open("lmStats.csv", "a+")
    tfp = open("lmRuntimes.csv", "a+")
    ofp = open(time_series_path + str(i) + "/lm.csv", "w")
    save_output_and_create_stats( inj,
                                  rs.stdout,
                                  sfp,
                                  ofp,
                                  ws,
                                  tfp,
                                  tb
                                )
  
    tb = time.time()
    #https://github.com/jMotif/SAX
    #EMMA: java -Xmx7g -jar EMMA.jar filename motifsize motifrange paasize alphabetsize snormalizationthreshold
    #java -Xmx7g -jar EMMA.jar timeSeriesTopLatentMotifBenchmark/time_series_0.csv 200 1.88 4 3 0.01
    rs = subprocess.Popen( [ "java",
                             "-Xmx7g",
                             "-jar",
                             "algorithms/EMMA.jar",
                             ts_abs_path,
                             str(ws),
                             str(r),
                             "6",
                             "4",
                             "0.01"
                           ],
                           stdout = subprocess.PIPE
                         )
    sfp = open("emmaStats.csv", "a+")
    tfp = open("emmaRuntimes.csv", "a+")
    ofp = open(time_series_path + str(i) + "/emma.csv", "w")
    save_output_and_create_stats( inj,
                                  rs.stdout,
                                  sfp,
                                  ofp,
                                  ws,
                                  tfp,
                                  tb
                                )
    #with 2*r
    rs = subprocess.Popen( [ "java",
                             "-Xmx7g",
                             "-jar",
                             "algorithms/EMMA.jar",
                             ts_abs_path,
                             str(ws),
                             str(2*r),
                             "6",
                             "4",
                             "0.01"
                           ],
                           stdout = subprocess.PIPE
                         )
    sfp = open("emma2rStats.csv", "a+")
    ofp = open(time_series_path + str(i) + "/emma_2r.csv", "w")
    save_output_and_create_stats( inj,
                                  rs.stdout,
                                  sfp,
                                  ofp,
                                  ws
                                )
  
    tb = time.time()
    #https://github.com/GrammarViz2/grammarviz2_src
    #GrammarViz3: java -Xmx7g -jar GrammarViz3.jar --data_in timeSeriesTopLatentMotifBenchmark/time_series_0.csv --window_size 200 --alphabet_size 4 --word_size 6 --strategy EXACT --threshold 0.01
    #java -Xmx7g -jar GrammarViz3.jar --data_in timeSeriesTopLatentMotifBenchmark/time_series_0.csv --window_size 200 --alphabet_size 4 --word_size 6 --strategy EXACT --threshold 0.01
    rs = subprocess.Popen( [ "java",
                             "-Xmx7g",
                             "-jar",
                             "algorithms/GrammarViz3.jar",
                             "--data_in",
                             ts_abs_path,
                             "--window_size",
                             str(ws),
                             "--alphabet_size",
                             "4",
                             "--word_size",
                             "6",
                             "--strategy",
                             "EXACT",
                             "--threshold",
                             "0.01"
                           ],
                           stdout = subprocess.PIPE
                         )
    sfp = open("gvStats.csv", "a+")
    tfp = open("gvRuntimes.csv", "a+")
    ofp = open(time_series_path + str(i) + "/gv.csv", "w")
    save_output_and_create_stats( inj,
                                  rs.stdout,
                                  sfp,
                                  ofp,
                                  ws,
                                  tfp,
                                  tb
                                )
  
    tb = time.time()
    #ScanMK
    rs = subprocess.Popen( [ "java",
                             "-Xmx7g",
                             "-jar",
                             "algorithms/ScanMK.jar",
                             "d=" + ts_abs_path,
                             "w=" + str(ws),
                             "k=1",
                             "r=" + str(r)
                           ],
                           stdout = subprocess.PIPE
                         )
    sfp = open("scanmkStats.csv", "a+")
    tfp = open("scanmkRuntimes.csv", "a+")
    ofp = open(time_series_path + str(i) + "/scanmk.csv", "w")
    save_output_and_create_stats( inj,
                                  rs.stdout,
                                  sfp,
                                  ofp,
                                  ws,
                                  tfp,
                                  tb
                                )
    
    tb = time.time()
    #ClusterMK
    rs = subprocess.Popen( [ "java",
                             "-Xmx7g",
                             "-jar",
                             "algorithms/ClusterMK.jar",
                             "d=" + ts_abs_path,
                             "w=" + str(ws),
                             "k=1",
                             "r=" + str(r)
                           ],
                           stdout = subprocess.PIPE
                         )
    sfp = open("clustermkStats.csv", "a+")
    tfp = open("clustermkRuntimes.csv", "a+")
    ofp = open(time_series_path + str(i) + "/clustermk.csv", "w")
    save_output_and_create_stats( inj,
                                  rs.stdout,
                                  sfp,
                                  ofp,
                                  ws,
                                  tfp,
                                  tb
                                )
    
    tb = time.time()
    #SetFinder
    rs = subprocess.Popen( [ "java",
                             "-Xmx7g",
                             "-jar",
                             "algorithms/SetFinder.jar",
                             "d=" + ts_abs_path,
                             "w=" + str(ws),
                             "k=1",
                             "r=" + str(r)
                           ],
                           stdout = subprocess.PIPE
                         )
    sfp = open("setfinderStats.csv", "a+")
    tfp = open("setfinderRuntimes.csv", "a+")
    ofp = open(time_series_path + str(i) + "/setfinder.csv", "w")
    save_output_and_create_stats( inj,
                                  rs.stdout,
                                  sfp,
                                  ofp,
                                  ws,
                                  tfp,
                                  tb
                                )
  
    inj.clear()
  
  if os.path.isdir("results"):
    shutil.rmtree("results")
  os.mkdir("results")
  for path in glob.glob("*Stats.csv"):
    shutil.move(path, "results")
  for path in glob.glob("*Runtimes.csv"):
    shutil.move(path, "results")

  #plot stats
  subprocess.Popen( [ "Rscript",
                      "generatePlots.R"
                    ],
                    stdout = subprocess.PIPE
                  )
  
  for path in glob.glob("*Plot.pdf"):
    shutil.move("*Plot.pdf", "results")
  
main()
