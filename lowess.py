# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 13:22:43 2014

@author: s4493222
"""
"""
https://code.google.com/p/crmeng-pre-1/source/browse/branches/start/lowess.py?r=4
This module implements the Lowess function for nonparametric regression.

Functions:
lowess        Fit a smooth nonparametric regression curve to a scatterplot.

For more information, see

William S. Cleveland: "Robust locally weighted regression and smoothing
scatterplots", Journal of the American Statistical Association, December 1979,
volume 74, number 368, pp. 829-836.

William S. Cleveland and Susan J. Devlin: "Locally weighted regression: An
approach to regression analysis by local fitting", Journal of the American
Statistical Association, September 1988, volume 83, number 403, pp. 596-610.
"""

import numpy
from numpy import median
import os
import string
SVNRevision="$Revision"
def lowess(x, y, f=2./3., iter=3):
    """lowess(x, y, f=2./3., iter=3) -> yest

Lowess smoother: Robust locally weighted regression.
The lowess function fits a nonparametric regression curve to a scatterplot.
The arrays x and y contain an equal number of elements; each pair
(x[i], y[i]) defines a data point in the scatterplot. The function returns
the estimated (smooth) values of y.

The smoothing span is given by f. A larger value for f will result in a
smoother curve. The number of robustifying iterations is given by iter. The
function will run faster with a smaller number of iterations."""
    n = len(x)
    r = int(numpy.ceil(f*n))
    #x = numpy.array(x)
    h = [numpy.sort(numpy.abs(x-x[i]))[r] for i in range(n)]
    w = numpy.clip(numpy.abs(([x]-numpy.transpose([x]))/h),0.0,1.0)
    w = 1-w*w*w
    w = w*w*w
    yest = numpy.zeros(n)
    delta = numpy.ones(n)
    for iteration in range(iter):
        for i in range(n):
            weights = delta * w[:,i]
            b = numpy.array([sum(weights*y), sum(weights*y*x)])
            A = numpy.array([[sum(weights), sum(weights*x)],
                             [sum(weights*x), sum(weights*x*x)]])
            beta = numpy.linalg.solve(A,b)
            yest[i] = beta[0] + beta[1]*x[i]
        residuals = y-yest
        s = numpy.median(abs(residuals))
        delta = numpy.clip(residuals/(6*s),-1,1)
        delta = 1-delta*delta
        delta = delta*delta
    return yest
    
class R:
  def __init__(self,Rpath='C:\\Program Files\\R\\R-3.0.2\\bin\\Rscript.exe', tmp="R_call_tmp.txt"):
    self.__Rpath='"'+Rpath + '" '
    self.__tmpfile=tmp
    self.__resultfn="return_"+tmp
    
  def ascsv(self, x, y, fn, Headings = None):
    with open(fn, 'w') as f:
      if Headings == None:
        f.write('"X", "Y"\n')
      else:
        head = ""
        for h in Headings:
          head.append('"')
          head.append(h)
          head.append('"')
        f.write( head+"\n")
      for i in range(min(len(x), len(y))):
        f.write("%s, %s\n" % (str(x[i]), str(y[i])))
    return fn
        
  def call(self, command):
    #print self.__Rpath + " " +command
    try:
      return os.system(self.__Rpath + command)
    except Exception as e:
      print "command", command
      raise e
  
  def lowess(self, inx, iny, f=2./3., iter=3, script=('C:\\Users\\s4493222\\Documents\\abrupt\\src\\projects\\sea_level_rise\\dolowess.R')):
    fn=self.ascsv(inx,iny, self.__tmpfile)
    command = script + " " + self.__tmpfile+" dummy dummy "+self.__resultfn+" " + str(f) 
    
    result = self.call(command) 
    x = []
    y=[]
    d1=[]
    d2=[]
    if result == 0:
      #valid command
      with open(self.__resultfn, "r") as res:
        for line in res.readlines()[1:]:
          (s, x1,y1, d0, dd1, dd2)=line.split()
          x.append(float(x1))
          y.append(float(y1))
          d1.append(float(dd1))
          d2.append(float(dd2))
          
      
    return numpy.array(inx), numpy.interp(inx, x, y),numpy.interp(inx, x, d1),numpy.interp(inx, x, d2)
    
  def EMD(self, inx, iny, detrend=True, boundary="wave", stop="type2"):
    return [[]]
    
  def strucchange(self, inx, iny, script=('C:\\Users\\s4493222\\Documents\\abrupt\\strucchange\\dostruct.R')):
    fn=self.ascsv(inx,iny, self.__tmpfile)
    command = script + " " + self.__tmpfile+" X Y "+self.__resultfn+" " + str(fn) 
    #print command
    result = self.call(command) 
    #print result
    if result == 0:
      #valid command
      with open(self.__resultfn, "r") as res:
        state = 0
        for line in res.readlines()[1:]:
          #print line
          if state == 0:
            pos = line.find("-segment partition")
            if pos > -1: #then the digits just before are the number of partitions
              pos2=line[:pos].rfind(" ") #the last space
              N = int(line[pos2:pos])
              #print "N",N
              state = 1
          elif state == 1:
            if line.find("% breakpoints") > -1:
              state = 2
              lread = 0
              result = []
          elif state == 2:
            #read N-1lines
            if lread < N-1:
             # print "Split",line
              (n, lo, brk, hi) = line.split()
              lread +=1
              result.append((int(lo),int(brk),int(hi)))
            else:
              state = 3
    return N, result
      
  def strucchange2way(self, inx, iny, script=('C:\\Users\\s4493222\\Documents\\abrupt\\strucchange\\dostructXY.R')):
    return self.strucchange(inx, iny, script=script)

  def changepoint(self, inx, iny, script=('C:\\Users\\s4493222\\Documents\\abrupt\\SpatialBPs\\dochangepoint.R')):
    fn=self.ascsv(inx,iny, self.__tmpfile)
    command = script + " " + self.__tmpfile+" X Y "+self.__resultfn+" " + str(fn) 
    #print command
    result = self.call(command)
    N = 0
    cpts = []
    #print result
    if result == 0:
      #valid command
      state = 0
      with open(self.__resultfn, "r") as res:
        lines = res.readlines()
        i = 0
        while state <2 and i < len(lines):          
          line = lines[i]
          if state == 0 and line.find('Slot "cpts":') > -1:
            state = 1
          if state == 1 and line.find("[1] ") > -1:
            cpts=line[4:].split()
            state = 2
          i += 1
      N = len(cpts)
    return N, [int(cp) for cp in cpts]
          
          
          
    
if __name__ == "__main__":
    import sys
    import window40test
    import numpy as np
    fn="C:\\Users\\s4493222\\Documents\\ReferenceData\\Rogers_Analysis_17Jul2014\\Data 4 Jim\\Global temps.csv"
    data=np.genfromtxt(fn,delimiter=",",names=True,filling_values =np.NaN) #fill with NaNs so must use their bounds
    print R().strucchange(data["Year"], data["CSIRO35A2r1"])
    print R().strucchange2way(data["CSIRO35A2r1"], data["Year"])
    print R().strucchange2way(data["Year"], data["CSIRO35A2r1"])
    print R().changepoint(data["Year"][21:], data["CSIRO35A2r1"][21:])
   # sys.exit()
#    all85=window40test.dataDict("C:\Users\\s4493222\\Documents\\abrupt\\CMIP5_breakpoints\\historical_rcp45qccceGW.txt")
#    model=all85.models()[0]
#    ys=all85[model]
#    xs=all85.years()[model]
#    #print xs, ys
#    x,y,d1,d2=R().lowess(xs, ys)
#    #print x,y,d1, d2
#    
#    #sys.exit()
#    
#    print "LOWESS",lowess(np.array(xs), np.array(ys),iter=1)
#    print R().call('""C:\\Users\\s4493222\\Documents\\abrupt\\src\\projects\\sea_level_rise\\dolowess.R"" ')
