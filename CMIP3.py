# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 10:43:37 2014

@author: s4493222
"""

#CMIP3.py reads in CMIP3 data for processing
import numpy as np

'''
File format looks as below
# Area weighted global warming (Deg K)
# Creation_date:  2007-05-13T14:56:32+10
# function: generateTasGW (ccir_lib.tcl: $Rev: 286 $)
# Source_file: '/cs/datastore/u/csdar/csar4/dawn_2100/SRESA1B_tas/pcmdi.ipcc4.csiro_mk3_5.sresa1b.run1.monthly.tas_A1_1871-2100.nc'
1871  287.812
'''
SVNRevision="$Revision: 286 $"
class CMIP3Exception(Exception):
  def __init__(self, msg):
    self.__msg="CMIP3Exeption:"+str(msg)
  def __str__(self):
    return str(self.__msg)
    
class CMIP3gw(object):
  def __init__(self, filename, sep=None):
    self.__Years=[]
    self.__Temperature=[]
    self.__metadata=[]
    self.__sep=sep
    self.__appendList__(filename)

  def __appendList__(self, filename):
    with open(filename, 'r') as f:
      lines = f.readlines()
      for line in lines:
        if line[0] == "#":
          self.__metadata.append(line)
        else:
          y,v=line.split(self.__sep)
          self.__Years.append(float(y))
          self.__Temperature.append(float(v))
    self.__Years=np.array(self.__Years, dtype=np.float32)
    self.__Temperature=np.array(self.__Temperature, dtype=np.float32)
    
  def appendFile(self, filename):    
    self.__Years=self.__Years.tolist()
    self.__Temperature=self.__Temperature.tolist()
    self.__appendList__(filename)
      
  TempConversion=['Raw', 'Kelvin', 'Celcius', 'Anomaly']
  
  
  def Years(self):
    return self.__Years
    
  def Temperature(self):
    return self.__Temperature
    
  def count(self):
    return len(self.__Years)
    
  def Warming(self, mode='Anomaly', baseincludes=[1975, 2004]):
    if mode == 'Anomaly':
      lo=self.__Years.tolist().index(baseincludes[0])
      hi=self.__Years.tolist().index(baseincludes[1])+1
      mean=self.__Temperature[lo:hi].mean()
      return self.__Temperature-mean
    elif mode== 'Kelvin' and self.__Temperature[0] < 200:
      return self.__Temperatures + 273.15
    elif mode=='Celcius' and self.__Temperature[0] > 200:
      return self.__Temperatures - 273.15
    elif not mode in self.TempConversion:
      raise CMIP3Exception('Conversion mode must be in '+str(self.TempConversion))
    else:
      return self.__Temperature

if __name__ == "__main__":
  fn='C:\\\Users\\s4493222\\Documents\\ReferenceData\\CSIRO_Model_Set\\gwFiles\\csiro_mk3_0.sresa1b.run1.monthly.tas_A1_1871-2100_gw.txt'
  fn='C:\\\Users\\s4493222\\Documents\\ReferenceData\\CSIRO_Model_Set\\gwFiles\\csiro_mk3_5.sresa1b.run1.monthly.tas_A1_1871-2100_gw.txt'

  c3=CMIP3gw(fn)
  for i in range(c3.count()):
    print c3.Years()[i], c3.Temperature()[i], c3.Warming()[i]
  
  import ConvergentBreaks as convergent_breaks
  import random
  ys =  c3.Warming()
  xs = np.array([random.random() for y in ys])
  Years= c3.Years()
  convergent_breaks.TraceFile=""
  cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', shallow=False, mode="control", guide="Stability",screenpr=0.05, pr=0.01, trace=True)
  print cb