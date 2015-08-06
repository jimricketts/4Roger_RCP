# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 12:02:29 2014

@author: s4493222
"""

#CMIP5.py reads in CMIP5 data for processing
import numpy as np
from CMIP3 import CMIP3gw
SVNRevision="$Revision: 307 $"
'''
File format looks as below (no headers, and comma separators)
1871,  287.812
'''
SVNRevision="$Revision: 307 $"
class CMIP5Exception(Exception):
  def __init__(self, msg):
    self.__msg="CMIP5Exeption:"+str(msg)
  def __str__(self):
    return str(self.__msg)
    
class CMIP5gw(CMIP3gw):
  def __init__(self, filenameList):
    super(CMIP5gw,self).__init__(filenameList[0], sep=',')
    for filename in filenameList[1:]:
      super(CMIP5gw,self).appendFile(filename)


if __name__ == "__main__":
  fn='C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_MIROC5_historical_r1i1p1_185001-201212.GW'
  c3=CMIP5gw([fn])
  for i in range(c3.count()):
    print c3.Years()[i], c3.Temperature()[i], c3.Warming()[i]
  
  import convergent_breaks
  import random
  ys =  c3.Warming()
  xs = np.array([random.random() for y in ys])
  Years= c3.Years()
  cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=0.05, pr=0.01, trace=False)
  print cb
