#! /usr/bin/python2.7

__author__='Diogo Vila Vicosa'
__version__='2016.05.03'

import argparse
import numpy as np

def ParseArguments():
  parser = argparse.ArgumentParser()
  parser.add_argument('-f',   action='store', dest='f',   help='COM dists file')
  parser.add_argument('-win', action='store', dest='win', help='Window size (in nm)')
  parser.add_argument('-low', action='store', dest='low', help='Lower distance (in nm)')
  parser.add_argument('-lar', action='store', dest='lar', help='Larger distance (in nm)')
  return parser.parse_args()

def closest_ndx(array,val):
  ndx=0
  ref=abs(array[0]-val)
  for i in range(len(array)):
    if abs(array[i]-val) < ref:
      ndx=i
      ref=abs(array[i]-val)
  return ndx

if __name__ == '__main__':
  inp = ParseArguments()
  low=float(inp.low)
  lar=float(inp.lar)
  win=float(inp.win)
  windows = np.arange(low,lar+win,win)
  #Read COM distances
  time = []
  dist = []
  with open(inp.f) as F:
    for l in F:
      parts = l.split()
      time.append(float(parts[0]))
# if z distance, field 3
      dist.append(float(parts[1]))
  #Select conformations
  ndxs = [closest_ndx(dist,i) for i in windows]
  for i in ndxs:
    print int(time[i]),dist[i]
