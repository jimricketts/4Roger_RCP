# -*- coding: utf-8 -*-
"""
Created on Mon Sep 08 15:22:58 2014

@author: s4493222
"""

import numpy as np
import bivariate_multi as bivariate
import scipy.stats as stats
import STARS
import ConvergentBreaks as convergent_breaks
import random
import whitening
import datetime
import copy
SVNRevision="$Revision: 286 $"

def shuffle_stat(xs1, ys1, Years, pos, iterations=100):
  '''
  recompute Ti0 and breakpoint locations by reanalysing the values on each side of the breaks point
  xs are control
  ys are test
  pos is the location being tested
  iterations is the number of times to test  
  '''
  #MUST take copies of input because the shuffles have side effects. JHR 13/10/2014
  if len(xs1) == 0:
    print "oops"
  xs=copy.copy(xs1)
  ys=copy.copy(ys1)  
  index=range(len(xs))
  TiList = []
  TiPosIndex=[]
  ShiftList=[]
  pos = min(pos, len(index))
  for i in range(iterations):
    lowx=xs[:pos+1]
    lowy=ys[:pos+1]
    hix=xs[pos+1:]
    hiy=ys[pos+1:]
    np.random.shuffle(lowx)
    np.random.shuffle(lowy)
    np.random.shuffle(hix)
    np.random.shuffle(hiy)
    lowx = np.append(lowx, hix)
    lowy = np.append(lowy, hiy)
    bv=bivariate.bivariate(lowy, lowx, anomalise=False, constantsxy=True)
    try:
      TiList.append(bv.maxTi())
      TiPosIndex.append(Years[bv.maxIndexTi()])
      ShiftList.append(bv.stepChange())
    except Exception as e:
      print str(e)
      raise
  return TiPosIndex, TiList,ShiftList
  

def shuffle_cut(xs, ys, Years, pos, iterations=100):
  '''
  recompute Ti0 and breakpoint locations by reanalysing the values on each side of the breaks point
  xs are control
  ys are test
  pos is the location being tested
  iterations is the number of times to test  
  '''
  
  TiPosIndex, TiList,ShiftList=shuffle_stat(xs, ys, Years, pos, iterations=iterations)
  bins=np.bincount(TiPosIndex)
  mode=bins.argmax()
  return stats.norm.fit(TiPosIndex), stats.norm.fit(TiList),stats.norm.fit(ShiftList), float(bins[mode])/iterations, mode, 
  
def shuffle_mode(xs, ys, Years, pos, iterations=100):
  '''
  recompute Ti0 and breakpoint locations by reanalysing the values on each side of the breaks point
  xs are control
  ys are test
  pos is the location being tested
  iterations is the number of times to test  
  '''
  
  TiPosIndex, TiList,ShiftList=shuffle_stat(xs, ys, Years, pos, iterations=iterations)
  bins=np.bincount(TiPosIndex)
  mode=bins.argmax()
  return bins[mode], mode

def jack_knife(xs, ys, Years, pos, iterations=100):
  ylist=[]
  mlist=[]
  for i in range(1,len(Years)-1):
    count, year=shuffle_mode(xs, ys, Years, i, iterations)
    ylist.append(year)
    mlist.append(count)
  bins=np.bincount(ylist)
  mode=bins.argmax()
  return float(bins[mode])/iterations, mode
  
  
  
if __name__ == "__main__":
  prewhiten=False
  outf=open("SHUFFLE.txt", "w")
  print >>outf, datetime.datetime.now()
  data=np.genfromtxt(STARS.fn,delimiter=",",names=True,filling_values =np.NaN) #fill with NaNs so must use their bounds
  for d in ["A4"]:#data.dtype.fields.keys():
    ys = data[d]
    if prewhiten:
      ys=whitening.prewhiten(ys,STARS.AlphaEst(ys, 15, option="optIPN4", returnmsgs=False))
    Years=data['Year']
    xs= np.array([random.random() for i in ys])
    for i in range(len(ys)):
      print >>outf, Years[i],xs[i], ys[i]
    cb=convergent_breaks.convergentBreaks(ys, xs, Years, d, screenpr=0.01, pr=0.01, trace=False)
    breaks=[0]
    for b in cb[2]:
      #print b
      breaks.extend([list(Years[:]).index(round(b[0][0])) +1])
    breaks.append(len(Years))
    
    #print "BREAKS",breaks
    for b in range(len(breaks))[1:-1]:
      sc=shuffle_cut(xs[breaks[b-1]:breaks[b+1]], ys[breaks[b-1]:breaks[b+1]], Years[breaks[b-1]:breaks[b+1]], breaks[b]-breaks[b-1])
      #jk=jack_knife(xs[breaks[b-1]:breaks[b+1]], ys[breaks[b-1]:breaks[b+1]], Years[breaks[b-1]:breaks[b+1]], breaks[b]-breaks[b-1])
      print >>outf, "CUT",d, b, breaks[b], cb[2][b-1][0], sc
      
  outf.close()
  #an experiment which simply quantifies the effect of a break at half  trend of 0.5, 1.0, 1.5 2.0, 2.5 3.0 Std 
  std=stats.tstd(xs)
  length=len(xs)
  factor = 0.5 
  with open("TREND_Break_%f.txt" % (factor,),"w") as outf:
    for k in range(10000):
      ys=np.array([random.random() + factor*std*i/length for i in range(length)])
      for i in range(int(length / 2)):
        ys[i] -= factor*std*i/length
      cb=convergent_breaks.convergentBreaks(ys, xs, Years, str(factor), screenpr=0.01, pr=0.01, trace=False)
      breaks=[0]
      for b in cb[2]:
        #print b
        breaks.extend([list(Years[:]).index(round(b[0][0])) +1])
      breaks.append(len(Years))
      for b in range(len(breaks))[1:-1]:
        sc=shuffle_cut(xs[breaks[b-1]:breaks[b+1]], ys[breaks[b-1]:breaks[b+1]], Years[breaks[b-1]:breaks[b+1]], breaks[b]-breaks[b-1])
        #jk=jack_knife(xs[breaks[b-1]:breaks[b+1]], ys[breaks[b-1]:breaks[b+1]], Years[breaks[b-1]:breaks[b+1]], breaks[b]-breaks[b-1])
        print >>outf, k, "CUT",d, b, breaks[b], cb[2][b-1][0:2], sc
        print k, "CUT",d, b, breaks[b], cb[2][b-1][0:2], sc
  outf.close()
  #an experiment which simply quantifies the effect of a trend of 0.5, 1.0, 1.5 2.0, 2.5 3.0 Std 
#  std=stats.tstd(xs)
#  length=len(xs)
#  factor = 0.5 * 4
#  with open("TREND_%f.txt" % (factor,),"w") as outf:
#    for k in range(10000):
#      ys=np.array([random.random() + factor*std*i/length for i in range(length)])
#      cb=convergent_breaks.convergentBreaks(ys, xs, Years, str(factor), screenpr=0.01, pr=0.01, trace=False)
#      breaks=[0]
#      for b in cb[2]:
#        #print b
#        breaks.extend([list(Years[:]).index(round(b[0][0])) +1])
#      breaks.append(len(Years))
#      for b in range(len(breaks))[1:-1]:
#        sc=shuffle_cut(xs[breaks[b-1]:breaks[b+1]], ys[breaks[b-1]:breaks[b+1]], Years[breaks[b-1]:breaks[b+1]], breaks[b]-breaks[b-1])
#        #jk=jack_knife(xs[breaks[b-1]:breaks[b+1]], ys[breaks[b-1]:breaks[b+1]], Years[breaks[b-1]:breaks[b+1]], breaks[b]-breaks[b-1])
#        print >>outf, k, "CUT",d, b, breaks[b], cb[2][b-1][0:2], sc
#        print k, "CUT",d, b, breaks[b], cb[2][b-1][0:2], sc
#  outf.close()
        