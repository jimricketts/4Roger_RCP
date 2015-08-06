# -*- coding: utf-8 -*-
"""
Created on Fri Jan 16 16:28:47 2015

@author: s4493222
"""

import numpy as np
import regress
SVNRevision="$Revision: 385 $"

class brkrpt(object):
  def __init__(self, fn):
    self.filename=fn
    with open(fn, 'r') as f:
      lines=f.readlines()[1:]
    self.__data=[]
    state = 0
    for line in lines:
      if state == 0:
        if line[:6] == "BYEARS":
          state= 1
          self.__initial = eval("'"+line[6:-1]+"'")
          self.__data=np.array(self.__data,dtype=float)
          self.__trace=[]
        else:
          self.__data.append(np.array(line.split(), dtype=float))
      elif state ==1:
        if line[:9] =="Returning":
          state = 2
          pos1= line.find(' -> ')
          line2=line[9:pos1]
          line=line[pos1+4:]
          pos1=line.find('] [')
          s1=line[:pos1+1]
          s2=line[pos1+2:]
          self.__breaks=eval(s1)
          self.__stats =eval(s2[:-1])
          self.__recursive=eval(line2)
        else:
          self.__trace.append(line)
  def recursive(self):
    return self.__recursive
  def data(self):
    return self.__data
  def trace(self):
    return self.__trace
  def breaks(self):
    return self.__breaks
  def stats(self):
    return self.__stats 
  def years(self):
    return self.__data[:,3]
  def ys(self):
    return self.__data[:,1]
  def xs(self):
    return self.__data[:,2]
  def segments(self):
    tbreaks=self.breaks()
    yrs=list(self.years())
    los=[ yrs.index(b)+1 for b in tbreaks[:-1]]
    his=[ yrs.index(b)+1 for b in tbreaks[1:]]
    #his[-1]+=1 #because the upper bound is the last data year
    los[0]-=1
    yss=[ self.ys()[los[i]:his[i]] for i in range(len(los)) ]
    yrss=[ self.years()[los[i]:his[i]] for i in range(len(los)) ]
    xss=[ self.xs()[los[i]:his[i]] for i in range(len(los)) ]
    return yrss, yss, xss
    
  def setbreaks(self,breaklist):
    self.__breaks=np.array(breaklist)
    
if __name__ == "__main__":
  import os
  bp=brkrpt(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\4Roger_Nature_SVN_264\\HadCRUT.4.2.0.0.annual_ns_avg//HadCRUT.4.2.0.0.annual_ns_avg.txt_0.trace")    
  yearsegs, ysegs, xsegs = bp.segments()
  pre98Years=yearsegs[-2]
  pre98temps=ysegs[-2]
  post98Years=yearsegs[-1]
  post98temps=ysegs[-1]
  
  for i in range(len(pre98Years)):
    print "Pre", i, pre98Years[i], pre98temps[i]

  for i in range(len(post98Years)):
    print "Post", i+len(pre98Years), post98Years[i], post98temps[i]



  beta1,alpha1=regress.regress(pre98temps, pre98Years)

  beta2,alpha2=regress.regress(post98temps, post98Years)
  
  print beta1, alpha1, beta2, alpha2
  
  allYears=[y for y in pre98Years]
  allYears.extend(post98Years)
  print allYears

  alltemps=[t for t in pre98temps]
  alltemps.extend(post98temps)
  
  crossyhat1=-np.array([alpha1+beta1* y for y in allYears]) +alltemps
  crossyhat2=-np.array([alpha2+beta2* y for y in allYears]) +alltemps
  
  print "A",crossyhat1, np.sum(crossyhat1)
  print "B",crossyhat2, np.sum(crossyhat2)
  print "C", crossyhat2-crossyhat1, sum(crossyhat2-crossyhat1)
  
  beta3,alpha3=regress.regress(np.array(alltemps),np.array(allYears))
  crossyhat3=-np.array([alpha3+beta3* y for y in allYears]) +alltemps
  
  for i in range(len(crossyhat1)):
    print allYears[i], alltemps[i], crossyhat1[i], crossyhat2[i], allYears[i]*beta1+alpha1,allYears[i]*beta2+alpha2, crossyhat3[i],allYears[i]*beta3+alpha3
  
  import bivariate
  import random
  bv=bivariate.bivariate(crossyhat1,np.array([random.random() for c in crossyhat1]), anomalise=False, pr=0.01)
  bv2=bivariate.bivariate(crossyhat2,np.array([random.random() for c in crossyhat2]), anomalise=False, pr=0.01)
  bv3=bivariate.bivariate(crossyhat3,np.array([random.random() for c in crossyhat3]), anomalise=False, pr=0.01)
  bv4=bivariate.bivariate(alltemps,np.array([random.random() for c in alltemps]), anomalise=False, pr=0.01)

  