# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 17:12:14 2014

@author: s4493222
"""
SVNRevision="$Revision: 307 $"
#stability.py - code for examining the stability of breakpoints, especially near the start of finish of data
import bivariate
import tests17Aug
import numpy as np
import random

#def stability(testData, controlData, yearData):

data = np.genfromtxt(tests17Aug.fn,delimiter=",",names=True,filling_values =np.NaN)
ys=data["B4"]
l=len(ys)
std=np.std(ys)
brk=l/2+20
ys[brk:]+=std/2

seg=brk-20

Years=data["Year"]
xs=np.array([random.random() for y in Years])

for seg in range(brk+10):
  bv=bivariate.bivariate(ys[seg:l-34], xs[seg:l-34], anomalise=False, pr=0.01)
  print seg, Years[seg+1], Years[brk+1], bv.stepChange(), Years[bv.maxIndexTi()+1+seg], std,  bv.maxTi(),bivariate.Pr2(bv.maxTi(),99), bv.critical()
