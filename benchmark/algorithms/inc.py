#!/usr/bin/python

import sys
import os
import shutil

tsf = "time_series_"

# in: begin int, end int, number of folders
def main():
  b = int(sys.argv[1])
  e = int(sys.argv[2])
  d = int(sys.argv[3])
  if b != e:
    p = range(0, d)
    
    if b < e:
      p = reversed(p)

    for i in p:
      if os.path.isdir(tsf + str(b + i)):
        if os.path.isfile(tsf + str(b + i) + "/" + tsf + str(b + i) + ".csv"):
          shutil.move(tsf + str(b + i) + "/" + tsf + str(b + i) + ".csv", tsf + str(b + i) + "/" + tsf + str(e + i) + ".csv")
        if os.path.isfile(tsf + str(b + i) + "/" + tsf + "meta_" + str(b + i) + ".csv"):
          shutil.move(tsf + str(b + i) + "/" + tsf + "meta_" + str(b + i) + ".csv", tsf + str(b + i) + "/" + tsf + "meta_" + str(e + i) + ".csv")
        if os.path.isfile(tsf + str(b + i) + "/" + tsf + "plot_" + str(b + i) + ".plt"):
          shutil.move(tsf + str(b + i) + "/" + tsf + "plot_" + str(b + i) + ".plt", tsf + str(b + i) + "/" + tsf + "plot_" + str(e + i) + ".plt")
        shutil.move(tsf + str(b + i), tsf + str(e + i))

main() 
