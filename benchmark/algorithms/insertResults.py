#!/usr/bin/python

import sys
import os
import shutil

rf = "results/"

def merge_data(fn, l):
  fpr = open("./results/" + fn, "r")
  fpc = open("./" + fn, "r")
  fpn = open("./new" + fn, "a+")

  i = 0

  while i < l:
    fpn.write(fpr.readline())
    i = i + 1

  for line in fpc.readlines():
    fpn.write(line)
  
  for line in fpr.readlines():
    fpn.write(line)
  
  fpr.close()
  fpc.close()
  fpn.close()

# in: line number int
def main():
  l = int(sys.argv[1])
  
  for alg in ["clustermk", "emma", "gv", "lm", "mk", "mp", "scanmk", "setfinder"]:
    if os.path.isfile("./new" + alg + "Stats.csv"):
      os.remove("./new" + alg + "Stats.csv")
    merge_data(alg + "Stats.csv", l)
    shutil.move("./new" + alg + "Stats.csv", "./results/" + alg + "Stats.csv")
    
    if os.path.isfile("./new" + alg + "Runtimes.csv"):
      os.remove("./new" + alg + "Runtimes.csv")
    merge_data(alg + "Runtimes.csv", l)
    shutil.move("./new" + alg + "Runtimes.csv", "./results/" + alg + "Runtimes.csv")
  
  if os.path.isfile("./newemma2rStats.csv"):
    os.remove("./newemma2rStats.csv")
  merge_data("emma2rStats.csv", l)
  shutil.move("./newemma2rStats.csv", "./results/emma2rStats.csv")

main() 
