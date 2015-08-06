# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 19:19:18 2015

@author: s4493222
"""
SVNRevision="$Revision: 307 $"
#ICMP file reader
import numpy as np
class ICMP(object):
  def __init__(self, fn):
    with open(fn,"r") as f:
      lines = f.readlines()
    self.__filename = fn
    self.__header=lines[:3]
    self.__data=[]
    for line in lines[3:]:
      self.__data.append(np.array(line.split(),dtype=float))
    self.__data=np.array(self.__data)
    #print np.shape(self.__data)

  def data(self):
    return self.__data
    
  def monthly(self):
    return self.__data[:,1:]
    
  def annual_means(self):
    return np.mean(self.monthly(),axis=1)

  def monthly_means(self):
    return np.mean(self.monthly(),axis=0)

  def annual_sums(self):
    return np.sum(self.monthly(),axis=1)

  def monthly_sums(self):
    return np.sum(self.monthly(),axis=0)

  def years(self):
    return self.__data[:,0]
 
  def filename(self):
    return self.__filename
if __name__ == "__main__":
  icmp=ICMP('C:\\\Users\\s4493222\\Documents\\ReferenceData\\icmip5\\r1i1p3_GISS-E2-R_icmip5_tas_Amon_ens_piControl_0-360E_-90-90N_n_su_027.dat')

