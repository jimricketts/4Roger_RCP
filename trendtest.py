# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 15:27:36 2015

@author: s4493222
"""
import numpy as np
import bivariate_multi as bivariate
import random
SVNRevision="$Revision: 308 $"
#self, rawdata, xs, anomalise=True, step=1, averagesteps=False, critical=None, pr= None, window = 5, constantsxy=True):


ys=np.array([random.random() + i/100. for i in range(100)])

xs=np.array([random.random() for i in range(100)])

bv=bivariate.bivariate(ys, xs, pr=0.05, anomalise=False, window=1)

print "Trended",bv.maxIndexTi(), bv.maxTi()

ys2=ys
#ys2=np.array([random.random() + 0.1 + int(i >49) * 0.1 for i in range(100)])
#ys2=np.array([random.random() + 0.1 + int(i >49) * 0.1 + i/100. for i in range(100)])
#ys2=np.array([random.random() + 0.1 + int(i >49) * 0.1 + i/10. for i in range(100)])
ys2=np.array([random.random() + 0.1 + int(i >49) * 0.1 + (i/100.) * int(abs(i-49) <3) for i in range(100)])
#ys2=ys

bv=bivariate.bivariate(ys2, xs, pr=0.05, anomalise=False, window=1)

print bv.maxIndexTi(), bv.maxTi()

loc=bv.maxIndexTi()
for i in range(-5,5):
  bv=bivariate.bivariate(ys2[loc-15+i:loc+15+i], xs[loc-15+i:loc+15+i], pr=0.05, anomalise=False, window=1)
  print i,loc+i, loc-15+bv.maxIndexTi()+i, bv.maxIndexTi(), bv.maxTi()
  
