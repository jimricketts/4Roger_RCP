# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 13:41:14 2015

@author: s4493222
"""
import numpy as np
SVNRevision="$Revision: 384 $
#HADCRUT file reader
#http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/series_format.html

#==============================================================================
# These 'best estimate' series are computed as the medians of regional time series computed for each of the 100 ensemble member realisations. Time series are presented as temperature anomalies (deg C) relative to 1961-1990.
# 
# Quoted uncertainties are computed by integrating across the distribution described by the 100 ensemble members, together with additional measurement and sampling error and coverage uncertainty information.
# 
# The data files contain 12 columns:
# 
#     Column 1 is the date.
#     Column 2 is the median of the 100 ensemble member time series.
#     Columns 3 and 4 are the lower and upper bounds of the 95% confidence interval of bias uncertainty computed from the 100 member ensemble.
#     Columns 5 and 6 are the lower and upper bounds of the 95% confidence interval of measurement and sampling uncertainties around the ensemble median. These are the combination of fully uncorrelated measurement and sampling uncertainties and partially correlated uncertainties described by the HadCRUT4 error covariance matrices.
#     Columns 7 and 8 are the lower and upper bounds of the 95% confidence interval of coverage uncertainties around the ensemble median.
#     Columns 9 and 10 are the lower and upper bounds of the 95% confidence interval of the combination of measurement and sampling and bias uncertainties.
#     Columns 11 and 12 are the lower and upper bounds of the 95% confidence interval of the combined effects of all the uncertainties described in the HadCRUT4 error model (measurement and sampling, bias and coverage uncertainties).
#==============================================================================

class HADCRUT4(object):
  def __init__(self, filename, annual=True):
    with open(filename, 'r') as f:
      self.__annual=annual
      self.__filename = filename
      lines=f.readlines()
      self.__data = []
      for line in lines:
        words = line.split()
        if not annual:
          v=words[0].split('/')
          words[0] = str(float(v[0])+(float(v[1])-0.5)/12.0)
        self.__data.append(np.array(words, dtype=float))
      self.__data = np.array(self.__data)
    self.Seasons=["DJF","MAM","JJA","SON"]
    self.SeasonMonths=[[12,1,2],[3,4,5],[6,7,8],[9,10,11]]
    
  def __getitem__(self, month=None):
    if month==None or self.__annual:
      return self.annual()
    else:
      return self.monthly(month)
    
  def data(self):
    return self.__data
    
  def annual(self):
    if self.__annual:
      return self.__data[:,1]
    else:
      leng=len(self.__data)
      nyr=int( leng /12)
      leng = nyr * 12
      return (self.__data[:,1][:leng]).reshape(nyr,12).mean(axis=1)
      
  def monthly(self,month):
    if self.__annual:
      return None
    else:
      leng=len(self.__data)
      nyr=int( leng /12)
      leng = nyr * 12
      return (self.__data[:,1][:leng]).reshape(nyr,12)[:,month]
    
  def years(self):
    if self.__annual:
      return self.__data[:,0]
    else:
      leng=len(self.__data)
      nyr=int( leng /12)
      leng = nyr * 12
      return np.array(self.__data[:leng:12,0], dtype=int)
 
  def filename(self, Label=None):
    if Label == None:
      return self.__filename
    else:
      return self.__filename+"_"+Label
  def dates(self):
    return self.__data[:,0]    
  def seasonal(self, column, season, firstmonth=12):
    #compose means according to the rule that numerically greater months preceding lesser months are taken from the prior year
      #essentially by insisting on the seasonal start month being defined.
    first=0
    ym=self.dates() #dates
    firstdate = ym[0]+(firstmonth-1)/12.0
    while abs(firstdate - ym[first]) > 0.0625:
      first +=1
    last = int((len(ym)-first)/12.0) * 12 + first
    trialdata=self.data()[first:last][:,1]
    
    monthlist=[]
    if season in self.Seasons:
      monthlist=self.SeasonMonths[self.Seasons.index(season)]
    elif type(season) is list:
      monthlist=season
    #adjust monthlist for extracted first date and for 1s offset
    for i in range(len(monthlist)):
      monthlist[i] = (((monthlist[i])+11 - first) % 12)
    
    trialview=trialdata.view()
    trialdata = trialview.reshape((len(trialview) / 12,12))[:,monthlist]
    return np.mean(trialdata, axis=1), np.array(self.dates()[first:last]).reshape((len(trialview) / 12,12))[:,monthlist]
      
'''
 #Yr    jan    feb ...                                                                      annual
 #Yr    %earth coverage by month
 1851  0.770  0.308 -0.631  0.503  0.011  0.036  0.346  0.255  0.126  0.654 -0.131 -0.269  0.165
 1851      3      3      3      3      3      3      3      3      3      3      3      3
 1852  0.091  0.113 -0.088 -0.917  0.382  0.465  0.529  0.054  0.134 -0.299  0.131  1.819  0.201
 1852      3      3      3      3      3      3      3      3      3      3      3      3
 1853  0.657 -0.689 -1.383 -0.510 -0.244 -0.232  0.382  0.297 -0.364  0.058 -0.554 -1.281 -0.322
 1853      3      3      3      3      3      3      3      3      3      3      3      3
'''
from itertools import islice
class CRU(object):
  def __init__(self, filename):
    with open(filename, 'r') as f:
      self.__filename = filename
      lines=f.readlines()
      self.__dataraw = []
      self.__coverageraw=[]
      datalines = islice(lines, 0, len(lines), 2)
      infolines = islice(lines, 1, len(lines), 2)
      for line in datalines:
        words = line[:-1].split()
        self.__dataraw.append(np.array(words, dtype=float))
      for line in infolines:
        words = line[:-1].split()
        self.__coverageraw.append(np.array(words, dtype=float))
      
      self.__dataraw = np.array(self.__dataraw)
      self.__coverageraw = np.array(self.__coverageraw)
      
      if self.__coverageraw[-1][-1] == 0:
        #then the year is not complete
        self.__data = self.__dataraw.view()[:-1]
        self.__coverage = self.__coverageraw.view()[:-1]
      else:
        self.__data = np.view(self.__dataraw)
        self.__coverage = np.view(self.__coverageraw)
    self.Seasons=["DJF","MAM","JJA","SON"]
    self.SeasonMonths=[[12,1,2],[3,4,5],[6,7,8],[9,10,11]]
      
  def __getitem__(self, month=None):
    return self.monthly(month)
    
  def data(self):
    return self.__data
    
  def annual(self):
    return self.__data[:,-1]
      
  def monthly(self,month):
    if month==None:
      return self.__data[:,-1]
    else:
      return self.__data[:,1+month]
    
  def coverage(self,month):
    return self.__coverage[:,1+month]

  def years(self, all=False):
    if all:
      return self.__dataraw[:,0]
    else:
      return self.__data[:,0]
 
  def filename(self, Label=None):
    if Label == None:
      return self.__filename
    else:
      return self.__filename+"_"+Label

  def dates(self):
    yrs=self.years(all=True)
    monthvals=[(i+0.5)/12.0 for i in range(12)]
#    datevals = [(yr + mv)  for yr in yrs for mv in monthvals]
    datevals = np.array([yr + monthvals  for yr in yrs]).flatten()
    covval=self.__coverageraw[:,1:].flatten()
    first=0
    while covval[first] == 0.:
      first += 1
    last = len(covval)-1
    while covval[last] == 0.:
      last -= 1
    last+=1
    return np.array(datevals[first:last])    
    
  def seasonal(self, column, season, firstmonth=12):
    #compose means according to the rule that numerically greater months preceding lesser months are taken from the prior year
      #essentially by insisting on the seasonal start month being defined.
    first=0
    ym=self.dates() #dates
    firstdate = ym[0]+(firstmonth-1)/12.0
    while abs(firstdate - ym[first]) > 0.0625:
      first +=1
    last = int((len(ym)-first)/12.0) * 12 + first
    trialdata=self.__dataraw[:,1:-1].flatten()[first:last]
    
    monthlist=[]
    if season in self.Seasons:
      monthlist=self.SeasonMonths[self.Seasons.index(season)]
    elif type(season) is list:
      monthlist=season
    #adjust monthlist for extracted first date and for 1s offset
    for i in range(len(monthlist)):
      monthlist[i] = (((monthlist[i])+11 - first) % 12)
    
    trialview=trialdata.view()
    trialdata = trialview.reshape((len(trialview) / 12,12))[:,monthlist]
    return np.mean(trialdata, axis=1), np.array(ym[first:last]).reshape((len(trialview) / 12,12))[:,monthlist]

if __name__ == "__main__":
#  fn="C:\\Users\\s4493222\\Documents\\ReferenceData\\Hadley_Data\\18Mar2015_home\\HadSST.3.1.1.0_monthly_globe_ts.txt"
#  hc=HADCRUT4(fn,annual=False)
#  print hc.annual()
#  print hc.years()
#  for i in range(12): 
#    print i, hc[i]
#  print hc.seasonal("","SON") 
    
  fn="C:\\Users\\s4493222\\Documents\\ReferenceData\\Hadley_Data\\18Mar2015_home\\FromCRU\\CRUTEM4-gl.dat.txt"
  hc=CRU(fn)
  print hc.seasonal("","DJF")
  print hc.annual()
  print hc.years()
  #for i in range(12): print i, hc[i]
  
