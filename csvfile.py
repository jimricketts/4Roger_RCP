# -*- coding: utf-8 -*-
"""
Created on Mon Oct 06 16:03:07 2014

@author: s4493222
"""

fn="C:\\Users\\s4493222\\Documents\\ReferenceData\\Rogers_Analysis_17Jul2014\\\Data 4 Jim\\Global temps.csv"

SVNRevision="$Revision: 312 $"
import numpy as np

class CSVfile(object):
  def __init__(self, filename, delimiter=",", filling_values =np.NaN, timevar="Year", skip_header=0):
    self.__header=[]
    self.__filename=filename
    if skip_header >0:
      #cheat by opening the file twice
      with open(filename,'r') as f:
        self.__header=f.readlines()[:skip_header]
    self.__data=np.genfromtxt(filename,delimiter=delimiter,names=True,filling_values =filling_values,skip_header=skip_header) #fill with NaNs so must use their bounds
    self.__time=timevar
    self.__fields=[name for name in np.sort(self.__data.dtype.fields.keys())]
    self.__ranges={}
    for f in self.__fields:
      self.__ranges[f]=self.valid_range(f)
      
  def valid_range(self, name):
    lo=0
    hi=len(self.__data)
    if hi >0:
      while np.isnan(self.__data[name][hi-1]):
        hi -= 1
      while np.isnan(self.__data[name][lo]) and lo < hi:
        lo += 1
    if lo >= hi:
      lo = 0 
      hi = -1
    return (lo, hi)
  
  def __getitem__(self,key):
    return self.__data[key]
    
  def data4name(self, name):
    (lo,hi)=self.__ranges[name]
    return self.__data[name][lo:hi]
    
  def time4name(self, name):
    (lo,hi)=self.__ranges[name]
    return self.__data[self.__time][lo:hi]
  
  def fields(self):
    return self.__fields
  
  def header(self):
    return self.__header
   
  def filename(self):
    return self.__filename
if __name__ == "__main__":
  a = CSVfile(fn)
  print a.fields()
  print a.data4name('GFDL2A1r1')
    
    

