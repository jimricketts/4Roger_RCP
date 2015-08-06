# -*- coding: utf-8 -*-
"""
Created on Tue Oct 07 16:35:42 2014

@author: s4493222

Code to read the following format NOAA ascii
file name convention for areal average (aravg) time series:
ann=annual average
mon=monthly average
land_ocean=merged land-ocean surface temperature
land=land surface temperature
ocean=ocean surface temperature
latitudes=southern and northern limits of areal average
v=version number
yyyymm=date for the latest data

Annual data (aravg.ann.*) :
1st column = year
2nd column = anomaly of temperature (K)
3rd column = total error variance (K**2)
4th column = high-frequency error variance (K**2)
5th column = low-frequency error variance (K**2)
6th column = bias error variance (K**2)

Monthly data (aravg.mon.*) :
1st column = year
2nd column = month
3rd column = anomaly of temperature (K)
4th column = total error variance (K**2)
5th column = high-frequency error variance (K**2)
6th column = low-frequency error variance (K**2)
7th column = bias error variance (K**2)
8th column = diagnostic variable
9th column = diagnostic variable
10th column= diagnostic variable

NOTE: anomalies are based on the climatology from 1971 to 2000
"""

import datetime
import numpy as np

def ExceptionInfo():
    import linecache
    import sys
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    return 'EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj)

  
SVNRevision="$Revision: 383 $"
class NOAAascii(object):
  def __init__(self,fn,annual=False):
    self.__datadic={}
    self.__annual=annual
    try:
      if annual:
        self.__header=['year', 'anomt', 'vart', 'hfvar', 'lfvar', 'biasvar']
        with open(fn,'r') as f:
          lines=f.readlines()
        for line in lines[:]:
          year, anomt, vart, hfvar, lfvar, biasvar= line.split()
          year=int(year) 
          anomt=float(anomt)
          vart=float(vart)
          hfvar=float(hfvar)
          lfvar=float(lfvar)
          biasvar=float(biasvar)
          date=datetime.date(year,6,30)
          if not year in self.__datadic:
            self.__datadic[year] = {}
            self.__datadic[year][0]=[]
          self.__datadic[year][0].append((year, anomt, vart, hfvar, lfvar, biasvar))
      else:  
        self.__header=['year', 'month', 'anomt', 'vart', 'hfvar', 'lfvar', 'biasvar', 'diag1', 'diag2', 'diag3']
        with open(fn,'r') as f:
          lines=f.readlines()
        for line in lines[:]:
          year, month, anomt, vart, hfvar, lfvar, biasvar, diag1, diag2, diag3= line.split()
          year=int(year) 
          month = int(month) 
          anomt=float(anomt)
          vart=float(vart)
          hfvar=float(hfvar)
          lfvar=float(lfvar)
          biasvar=float(biasvar)
          diag1=float(diag1)
          diag2=float(diag2)
          diag3=float(diag3)
          date=datetime.date(year,month,15)
          if not year in self.__datadic:
            self.__datadic[year] = {}
          if not month in self.__datadic[year]:
            self.__datadic[year][month]=[]
          self.__datadic[year][month].append((year, month, anomt, vart, hfvar, lfvar, biasvar, diag1, diag2, diag3))
    except Exception as e:
      ei=ExceptionInfo()
      print ei
      print line
      raise Exception(str(ei))
  def monthly(self, variable='anomt', months=None, years=None, raw=False):
    result=[]
    count=0
    col = self.__header.index(variable)
    #colate variable by year and month
    for year in np.sort(self.__datadic.keys()):
      for month in np.sort(self.__datadic[year].keys()):
        
        count = len(self.__datadic[year][month])
        if variable in ['rain', 'evap']:
          value = np.sum([entry[col] for entry in self.__datadic[year][month]])
        else:
          value = np.mean([entry[col] for entry in self.__datadic[year][month]])
        if (years == None or year in years) and (months == None or month in months):
          result.append((year, month, value, count))
    if raw:
      return [r[2] for r in result]
    else:
      return result
    
  def annually(self, variable='anomt', years=None, raw=False, entireYears=True):
    result=[]
    count=0
    col = self.__header.index(variable)
    if self.__annual:
      for year in np.sort(self.__datadic.keys()):
        value = np.sum([entry[col] for entry in self.__datadic[year][0]])
        #this should only count complete years JHT 27May2015 Ticket #3
        result.append((year, 0, value, 1))
      if raw:
        return [r[2] for r in result]
      else:
        return result
      
    else:
      #colate variable by year and month
      for year in np.sort(self.__datadic.keys()):
        value=0.0
        count=0
        for month in np.sort(self.__datadic[year].keys()):
          value += np.sum([entry[col] for entry in self.__datadic[year][month]])
          count += len(self.__datadic[year][month])
        #this should only count complete years JHT 27May2015 Ticket #3
        if ((entireYears and count ==12) or (not entireYears and count > 0)) and (years==None or year in years):
          if variable in ['rain', 'evap']:
            result.append((year, 0, value, count))
          else:
            result.append((year, 0, value/count, count))
      if raw:
        return [r[2] for r in result]
      else:
        return result

if __name__ == "__main__":
  fn="C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\ZonalLandSeaASCII_NCDC\\aravg.mon.land_ocean.90S.90N.v3.5.4.201408.asc"
  a=NOAAascii(fn)
  print a.annually()
  print a.annually(entireYears=False)
