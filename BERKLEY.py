# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 18:27:57 2015

@author: s4493222
"""
SVNRevision="$Revision: 307 $"
import numpy as np
class BERKLEY(object):
  def __init__(self, fn="C://Users//s4493222//Documents//ReferenceData//Berkley//Land_and_Ocean//Land_and_Ocean_summary.txt"):
    with open(fn,'r') as f:
      self.__filename=fn
      lines = f.readlines()
      self.__header=[]
      self.__data=[]
      i=0
      while len(lines[i]) == 0 or lines[i][0] in ["#","%"]:
        self.__header.append(lines[i])
        i+=1
        #print lines[i]
      for line in lines[i:]:
        if len(line)>3:self.__data.append(np.array(line.split(),dtype=float))
      self.__data=np.array(self.__data, dtype=float)
      
  def data(self):
    return self.__data
    
  def years(self):
    return self.__data[:,0] 
    
  def annual(self):
    return self.__data[:, 1]
 
  def filename(self):
    return self.__filename
    
if __name__ == "__main__":
  bd=BERKLEY()    
