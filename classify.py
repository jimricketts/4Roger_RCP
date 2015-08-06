# -*- coding: utf-8 -*-
"""
Created on Wed Mar 04 21:05:19 2015

@author: James
"""

#classify_breaks
import numpy as np;
import bivariate_multi as bivariate
import random
def classify(ys, years, Year, window=10, span=2):
  spins=range(-span, span+1)
  counts=np.zeros(np.shape(ys))
  pos=list(years).index(Year)
#  if #need to think about what if window - span < 0
  for i in spins:
    lo=max(0, pos-window+i)
    hi=min(len(ys)+1, pos+window+i)
    yslice=ys[lo:hi]
#    yrslice=years[lo:hi]
#    slpos=list(yrslice).index(Year)
    for j in range(int(100/(2*span+1))):
      bv=bivariate.bivariate(yslice, np.array([random.random() for y in yslice]), anomalise=False, pr=0.01) #(testdata,controldata, anomalise=False, pr=0.01)
      bvpos=bv.maxIndexTi()
      counts[lo+bvpos]+=1
  return np.max(counts)/np.sum(counts)
  
if __name__ == "__main__":
  ys=np.array([0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1])
  ys=ys+np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9])
  years=np.array([1900+i for i in range(len(ys))])
  bv1=bivariate.bivariate(ys, np.array([random.random() for y in ys]), anomalise=False, pr=0.01)     
  print classify(ys, years, years[bv1.maxIndexTi()])   