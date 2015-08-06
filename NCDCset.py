# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 14:15:05 2015

@author: s4493222
"""

import numpy as np
import os
import glob
import csvfile

SVNRevision="$Revision: 313 $"
class NCDCFileSet(object):
  def __init__(self, fpath, keylist=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
    files=glob.glob(fpath)
    self.__years=[]
    self.__data=[]
    self.__keys=keylist
    self.__path=fpath
    for fn in files:
      for k in keylist:
        if fn.find(k) >=0:
          csv=csvfile.CSVfile(fn, timevar="Year",skip_header=2)
          if self.__years == []:
            self.__years=csv.data4name('Year')
          self.__data.append(csv.data4name('Value'))
    #trim off excess data
    leng=len(self.__years)
    for i in range(len(self.__data)):
      leng=min(leng, len(self.__data[i]))
    for i in range(len(self.__data)):
      if leng < len(self.__data[i]):
        self.__data[i]=self.__data[i][:leng]
    self.__data = np.array(self.__data)
    self.__years= np.array(self.__years[:leng])
    
  def __getitem__(self, i):
    return self.__data[i]
    
  def years(self):
    return self.__years
    
  def monthly(self, month):
    return self.__data[month-1]
    
  def annual(self):
    return np.mean(self.__data,0)
  
  def name(self, i):
    return self.__keys[i]
  
  def filename(self, month=None):
    pos1=self.__path.find('*')
    if pos1<0:
      return self.__path
    elif month==None:
      return self.__path[:pos1]+"-"+self.__keys[0]+"-"+self.__keys[-1]+"-"+self.__path[pos1+1:]
    else:
      return self.__path[:pos1]+"-"+self.name(month)+"-"+self.__path[pos1+1:]
    
if __name__ == "__main__":
  ncd=NCDCFileSet('C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_LandOcean\\SH_LandOcean*1880-2015.csv')
  print ncd.years()
  print ncd.annual()
  for i in range(1,13):
    print ncd[i-1]
  print ncd.filename()
  for i in range(0,12):
    print ncd.filename(i)   
