# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 10:27:22 2014

@author: s4493222
"""

import os
import numpy as np
#import bivariate
#import recursetest
import random
#import lowess
import datetime
#import matplotlib.pylab as plt
#import regress
#import math
#import copy
#import scipy.stats as stats
import tests17Aug
#import tests15Aug
#import whitening
import STARS


SVNRevision="$Revision: 307 $"
def prewhiten(x, rho):
  result = np.array(x)
  result[1:] -= rho * result[:-1]
  return result
  
if __name__ == "__main__":
  fn=os.environ["HOMEPATH"]+"\\Documents\\ReferenceData\\Rogers_Analysis_17Jul2014\\GISSTEMPV3.csv"
  fn=os.environ["HOMEPATH"]+"\\Documents\\abrupt\\CMIP5_breakpoints\\Artyfishul_dartar4Jim (1).csv"

  #fn = tests17Aug.fn
  data=np.genfromtxt(fn,delimiter=",",names=True,filling_values =np.NaN) #fill with NaNs so must use their bounds
  Years=data["Year"][:]       ####<<<<<<<<<<<<<<<<<<<<<<<<
  
  RunName = "whitening"
  TraceFile = RunName + "_Trace.txt"
  tf=open(TraceFile,"w")
  print >>tf,"trace", RunName, datetime.datetime.now()
  tf.close()
  
  f = open(RunName+".txt","w")
  
  discoveryList={} 
  initialList = {}
  reps=1
  cutoff=0.5
  for i in range(reps):
    Rand = np.array([random.random() for y in Years])

    print >>f,datetime.datetime.now()
    for p in data.dtype.names:
      model = p
      if not model in discoveryList:
        discoveryList[model]={}
      if not model in initialList:
        initialList[model]={}
      print >>f, model, RunName
      print "MODEL", model, i
      
    
      lo,hi=tests17Aug.valid(data[model][:])   #####<<<<<<<<<<<<<<<<<<<<<<<<<<
      #if model != "Year": 
      if model != "Year":
        try:
          Dat=data[model][lo:hi]#-sigmoid.sig4(Years[lo:hi],data[model][lo:hi])
        except:
          print "Cannot estimate sigmoid"
          Dat=data[model][lo:hi]
        brks=tests17Aug.convergentBreaks(Dat, Rand[lo:hi], Years[lo:hi], model, mode="control")
        try:
          Dat=prewhiten(data[model][lo:hi],STARS.AlphaEst(data[model][lo:hi], 15, option="optIPN4", returnmsgs=False))#-sigmoid.sig4(Years[lo:hi],data[model][lo:hi])
        except:
          print "Cannot estimate sigmoid"
          Dat=data[model][lo:hi]
        wbrks=tests17Aug.convergentBreaks(Dat, Rand[lo:hi], Years[lo:hi], model, mode="control")
        print  >>f,"MODEL",model, "FOUND", brks
        br=tuple(brks[1]) #immutable so can use as key
        wbr=tuple(wbrks[1])
        tests17Aug.RunName=RunName
        tests17Aug.plotbreaksSummary(Years[lo:hi], data[model][lo:hi], wbr, title=model+" "+os.path.basename(fn), 
                      subtitle="Breaks"+str(br[1:-1])+"\nWhiter"+str(wbr[1:-1]), breakyears=wbr, breakyearsAve=br, offset=0)
    
  #print discoveryList
  #print collate(discoveryList, cutoff=cutoff, reps=reps)
  f.close()
#    
