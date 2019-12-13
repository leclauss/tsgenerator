#!/usr/bin/python

import sys
import numpy

def sim(sub0, sub1):
  #calculate means
  meanPos0 = numpy.mean(sub0);
  meanPos1 = numpy.mean(sub1);
  
  stdPos0 = numpy.std(sub0);
  stdPos1 = numpy.std(sub1);
  
  #z-normalize subsequences
  normPos0 = []
  normPos1 = []
  
  for i in range(len(sub0)):
    normPos0.append((sub0[i] - meanPos0) / stdPos0)
    normPos1.append((sub1[i] - meanPos1) / stdPos1)
  
  #calculate similarity
  rs = 0.0
  
  for i in range(len(sub0)):
    rs += (normPos0[i] - normPos1[i]) * (normPos0[i] - normPos1[i])

  return numpy.sqrt(rs)

def main():
  
  if (len(sys.argv) != 5):
    print("Wrong number of arguments")
    exit(0)
  
  tsPath = sys.argv[1]
  ws = int(sys.argv[2])
  pos0 = int(sys.argv[3])
  pos1 = int(sys.argv[4])
  
  fp = open(tsPath, "r")
  
  ts = []
  
  for line in fp:
    ts.append(float(line))
  
  sub0 = ts[pos0:pos0+ws]
  sub1 = ts[pos1:pos1+ws]
  
  rs = sim(sub0, sub1) 
  
  print(numpy.sqrt(rs))

main()
