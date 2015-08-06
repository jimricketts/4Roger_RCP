# -*- coding: utf-8 -*-
"""
Created on Thu Sep 04 12:40:23 2014

@author: s4493222
"""
SVNRevision="$Revision: 286 $"

#APSIMread
#code for reading an APSIM file in from SILO and producing daily, montly and annual data


'''
[weather.met.weather]
!station number = 082039
!station name =  RUTHERGLEN RESEARCH
latitude = -36.1047  (DECIMAL DEGREES)
longitude =  146.5094  (DECIMAL DEGREES)
tav = 14.64 (oC) ! Annual average ambient temperature. Based on 1 Jan 1957 to current.
amp = 16.37 (oC) ! Annual amplitude in mean monthly temperature. Based on 1 Jan 1957 to current.
!Data Extracted from Silo on 20140904" for APSIM
! *** Some early data in this file is only possible because of the Climarc project, Thanks Climarc ***!As evaporation is read at 9am, it has been shifted to day before
!ie The evaporation measured on 20 April is in row for 19 April
!The 6 digit code indicates the source of the 6 data columns
!0 actual observation, 1 actual observation composite station
!2 interpolated from daily observations
!3 interpolated from daily observations using anomaly interpolation method for CLIMARC data
!6 synthetic pan
!7 interpolated long term averages
!more detailed two digit codes are available in SILO's 'Standard' format files
!
!For further information see the documentation on the datadrill
!  http://www.longpaddock.qld.gov.au/silo
!
year  day radn  maxt   mint  rain  evap    vp   code
 ()   () (MJ/m^2) (oC)  (oC)  (mm)  (mm) (hPa)     ()
1889   1   26.0  36.0  18.5  92.5   8.7  19.0 333263
'''

import datetime
import numpy as np
  

class APSIMfile(object):
  def __init__(self,fn):
    self.__datadic={}
    with open(fn,'r') as f:
      lines=f.readlines()
      start, header, self.__station, self.__name, self.__lat, self.__lon = self.splitAPSIMHeader(lines)
      self.__header=header[2:]
    for line in lines[start:]:
      year, day, radn, maxt, mint, rain, evap, vp, code= line.split()
      year=int(year) 
      day = int(day) 
      radn=float(radn)
      maxt=float(maxt)
      mint=float(mint)
      rain=float(rain)
      evap=float(evap)
      vp=float(vp)
      date=datetime.date(year,1,1)+datetime.timedelta(day-1)
      month=date.month
      if not year in self.__datadic:
        self.__datadic[year] = {}
      if not month in self.__datadic[year]:
        self.__datadic[year][month]=[]
      self.__datadic[year][month].append((radn, maxt, mint, rain, evap, vp))
      
  def splitAPSIMHeader(self, lines):
    i = 0
    while  i < len(lines) and lines[i][:22] != ' ()   () (MJ/m^2) (oC)  (oC)  (mm)  (mm) (hPa)     ()'[:22]:
      station = ""
      name = ""
      if lines[i][:18] == '!station number = ':
        station=lines[i][18:]
      if lines[i][:16] == '!station name = ':
        name=lines[i][16:]
      if lines[i][:10] == 'latitude =':
        lat=float(lines[i][10:lines[i].index('(')])
      if lines[i][:11] == 'longitude =':
        lon=float(lines[i][11:lines[i].index('(')])
      i += 1
    header = lines[i-1].split()
    return i+1, header, station, name, lat, lon
    
  def monthly(self, variable, months=None, years=None, raw=False):
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
    
  def annually(self, variable, years=None, raw=False):
    result=[]
    count=0
    col = self.__header.index(variable)
    #colate variable by year and month
    for year in np.sort(self.__datadic.keys()):
      value=0.0
      count=0
      for month in np.sort(self.__datadic[year].keys()):
        value += np.sum([entry[col] for entry in self.__datadic[year][month]])
        count += len(self.__datadic[year][month])
      if count > 0 and (years==None or year in years):
        if variable in ['rain', 'evap']:
          result.append((year, 0, value, count))
        else:
          result.append((year, 0, value/count, count))
    if raw:
      return [r[2] for r in result]
    else:
      return result

if __name__ == "__main__":
  cntdata=APSIMfile('C:\\Users\\s4493222\\Documents\\ReferenceData\\SILO\\rutheDRILL.sim.txt')
  tstdata=APSIMfile('C:\\Users\\s4493222\\Documents\\ReferenceData\\SILO\\rutheRGLEN.sim.txt')
  prewhiten=False
  variable='mint'
#  print np.mean(data.annually('maxt',years=[2012], raw=True))
  import convergent_breaks
  import random
  #cntdata.annually('maxt', years=Years, raw=True)

  Years=np.array(range(1910,2014))
  Rand = np.array([random.random() for r in Years])
  ys=np.array(tstdata.annually(variable, years=Years, raw=True))
  xs=np.array(cntdata.annually(variable, years=Years, raw=True))
  #ys=(ys+np.array(tstdata.annually('maxt', years=Years, raw=True)))/2.0
  
  if prewhiten:
    import STARS
    import whitening
    ys=whitening.prewhiten(ys,STARS.AlphaEst(ys, 15, option="optIPN4", returnmsgs=False))
    xs=whitening.prewhiten(xs,STARS.AlphaEst(ys, 15, option="optIPN4", returnmsgs=False))
  #xs=Rand
  print convergent_breaks.convergentBreaks(ys, xs, Years, "Rutherglen")