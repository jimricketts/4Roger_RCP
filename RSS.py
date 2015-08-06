# -*- coding: utf-8 -*-
"""
Created on Tue May 05 15:08:19 2015

@author: s4493222
"""
import numpy as np
#RSS satellite data
SVNRevision="$Revision: 346 $"
RSSMasterDir="C:\\Users\\s4493222\\Documents\\ReferenceData\\RSS\\Downloaded_20150505\\"

RSSfnPrefix="RSS_TS_channel_"
RSSfnPostfix="_Land_And_Sea_v03_3.txt"

RSSBands=["Continental%20US", "Global","Northern%20Hemisphere","Northern%20Mid%20Latitudes","Northern%20Polar","Northern%20Hemisphere","Northern%20Mid%20Latitudes","Northern%20Polar","Tropics"]
RSSdataset=['TLT','TMT','TLS',
            'TTS',
            'TTT']


RSSFiles=[RSSMasterDir+RSSfnPrefix+rds+"_" + rdb+ RSSfnPostfix for rds in RSSdataset for rdb in RSSBands]
#RSS_TS_channel_TLS_Continental%20US_Land_And_Sea_v03_3.txt


class RSS(object):
  def __init__(self, level):
    self.__metadata=[]
    self.__data=[]
    self.__dates={}
    self.__lofull=-1
    self.__hifull=-1
    self.__headers=[]
    self.__filenames=[]
    self.__level=level
    for (band,fn) in [(rdb,RSSMasterDir+RSSfnPrefix+rds+"_" + rdb+ RSSfnPostfix) for rds in [level] for rdb in RSSBands]:
      lines=[]
      self.__headers.append(band)
      self.__filenames.append(fn)
      row=0
      with open(fn, 'r') as f:
        lines = f.readlines()[:]
        self.__metadata=lines[:5]
        lines=lines[5:]
      for line in lines:
        (year, month, value)=line.split()
        year = int(year)
        month = int(month)
        value = float(value)
        date = year + (month-0.5)/12
        if value > -99.8:  
          if self.__lofull==-1:
            if month == 1:
              self.__lofull=date
          if self.__lofull > -1 and month==12:
            self.__hifull = date
          if not date in self.__dates:  
            self.__dates[date]=(year, month)
            self.__data.append([])
          self.__data[row].append(value)
          row += 1
    self.__data=np.array(self.__data)  
    keys=sorted(list(self.__dates.keys()))
    #print keys
    self.__firstfulldate=self.__lofull
    self.__lastfulldate=self.__hifull
    
    self.__lofull=keys.index(self.__lofull)
    self.__hifull=keys.index(self.__hifull)+1
    #print keys[self.__lofull:self.__hifull]

  def headers(self):
    return self.__headers

  def data(self, column=None,alldata=False):
    if column == None:
      if alldata:
        return self.__data
      else:
        return self.__data[self.__lofull:self.__hifull]
    else:
      if type(column) is int:
        col=column
      else:
        #print type(column)
        col = self.__headers.index(column)
      #print self.__lofull,self.__hifull
      if alldata:
        return self.__data[:,col]
      else:
        return self.__data[self.__lofull:self.__hifull,col]
  
  def yearrange(self):
    return range(int(self.__firstfulldate), int(self.__lastfulldate)+1)
    
  def annual(self, column):
    return  np.mean([self.data(column)[self.__lofull+mo:self.__hifull:12] for mo in range(12)],axis=0)  
  
  def filename(self, column):
    if type(column) is int:
      col = RSSBands[column]
    else:
      col = column
    return RSSMasterDir+RSSfnPrefix+self.__level+"_" + col+ RSSfnPostfix 
    
if __name__ == "__main__":
  rss=RSS('TLT')
  print rss.headers()
  print rss.annual("Continental%20US")   
