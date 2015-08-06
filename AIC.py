# -*- coding: utf-8 -*-
"""
Created on Tue Aug 05 10:22:39 2014

@author: s4493222
"""

#The Akaike Information Criterion Class - specifically for computing composite linear models.
#this version uses regress.analysed_regress
import regress
import numpy as np
import math
SVNRevision="$Revision: 307 $"
class AICException(Exception):
  def __init__(self, msg):
    self.__msg="AICExeption:"+str(msg)
  def __str__(self):
    return str(self.__msg)
    
class AIC(object):
  '''
  performa an AIC between testdata and control series broken by years corresponding to the subset of the datayears on the breaklist
  The main reason for using this is to enable multiple evaluations by saving the intermediate results indexed by segment bound
  testdata, controldata, datayears, breaks
  as of 17sep2014 this will also evaluate an AIC for steps only models
  '''
  def __init__(self, testdata, controldata, datayears, breaks, corrected=True, breaksonly=False):
    object.__init__(self)
    if len(breaks) < 2:
      raise AICException("break list must include start and end points (len >= 2)")
    self.__breaksonly=breaksonly  
    self.__segments = {}
    lo = np.argwhere(datayears==breaks[0])[0][0]
    hi = np.argwhere(datayears==breaks[-1])[0][0] +1
    if self.__breaksonly:      
      mean=np.mean(testdata[lo:hi])
      resid=(testdata[lo:hi]-mean)
      self.__segments[(lo, hi)]=(np.sum(resid*resid), hi-lo)#/mean/(hi-lo-1)
    else:
      stats = regress.analysed_regress(testdata[lo:hi], controldata[lo:hi])
      self.__segments[(lo, hi)]=(stats["SSE"], stats["n"])
    self.__result = self.evaluate( testdata, controldata, datayears, breaks, corrected)
  
  def value(self):
    return self.__result
    
  def segments(self):
    return self.__segments
    
  def evaluate(self, testdata, controldata, datayears,breaks, corrected=True):
    if len(breaks) < 2:
      raise AICException(msg="break list is incomplete in AIC.evaluate")
    try:
      #print "evaluate", breaks, datayears
      lo = np.argwhere(datayears==breaks[0])[0][0]
      hi = np.argwhere(datayears==breaks[-1])[0][0] +1
      statslist = []
      for b in range(1,len(breaks) -1):
        mid=np.argwhere(datayears==breaks[b])[0][0] # + 1 #
        key=(lo, mid)
        if not (key in self.__segments):
          if mid - lo <= 3:
            self.__segments[key]=(10.0e10, mid - lo -1)
          elif self.__breaksonly:
            mean=np.mean(controldata[lo:mid])
            resid=(controldata[lo:mid]-mean)/(mid-lo)#/mean
            self.__segments[key]=(np.sum(resid*resid), mid-lo)            
          else:
            stats = regress.analysed_regress(testdata[lo:mid], controldata[lo:mid])
            self.__segments[key]=(stats["SSE"], stats["n"])
        statslist.append(self.__segments[key])
        lo = mid
      key=(lo, hi)
      if not (key in self.__segments):
        if hi - lo <= 3:
          self.__segments[key]=(10.0e10, hi - lo-1 )
        elif self.__breaksonly:
          mean=np.mean(controldata[lo:hi])
          resid=(controldata[lo:hi]-mean)/(hi-lo-1)#/mean
          self.__segments[key]=(np.sum(resid*resid), hi-lo)            
        else:
          stats = regress.analysed_regress(testdata[lo:hi], controldata[lo:hi])
          self.__segments[key]=(stats["SSE"], stats["n"])
      statslist.append(self.__segments[key])
      return self.AIC(statslist, corrected=corrected)  
    except Exception as e:
      raise AICException(msg="breaklist " +str(breaks)+" gives Exception "+str(e))
      
  def evaluateLinear(self, testdata, controldata, datayears,breaks, corrected=True):
    return self.evaluate(testdata, controldata, datayears, np.array([breaks[0], breaks[-1]]), corrected)
    
  def AIC(self, statslist, corrected=True):
    '''
    Compute AIC for a composite model of regression lines.
    #and in a published paper 
    #https://www.researchgate.net/post/What_is_the_AIC_formula
    #Alternate AIC = N*log(RSS/N) + 2k is used
    #and suggested correction is applied by default to converge on the BIC behaviour
    '''
    if statslist == []:
      return np.NaN
    else:
      rss = 0
      n = 0
      k = len(statslist)
      if self.__breaksonly:
        k += 1  #only the number of segments plus 1 for the test
      else:
        k = 3 * k -1 #regression params plus break points
      for stats in statslist:
        rss += stats[0]
        n += stats[1]   
      try:
  #      return k + n * math.log( 2.0 * math.pi* rss/(n - k)) + 1, rss, n, k
        if corrected:
          correction=2*k*(k+1)/(n-k-1)
        else:
          correction=0.0
        return 2*k + n * math.log(rss/n) +correction , rss, n, k
      except Exception as e:
        raise AICException(str(e)+str(statslist))
        pass

if __name__ == "__main__":
  
  import numpy as np 
  import random
  fn="C:\\Users\\s4493222\\Documents\\ReferenceData\\Rogers_Analysis_17Jul2014\\GISSTEMPV3.csv"

  data=np.genfromtxt(fn,delimiter=",",names=True,filling_values =np.NaN) #fill with NaNs so must use their bounds
  
  Years=data["Year"][:]       ####<<<<<<<<<<<<<<<<<<<<<<<<
  Rand=np.array([random.random() for y in Years])
  
  
  for p in ["44S24S"]:#data.dtype.names:
  #for p in data.dtype.names:
    model = p    
      
    lo=0
    hi=len(data)
    a=AIC(data[model][lo:hi], Rand[lo:hi], Years[lo:hi], np.array([1880.0, 1892.0, 1931.0, 1968.0, 1970.0, 1995.0, 2012.0]))
    print "Linear Rand",a.value()
    print "Linear Rand",a.evaluate(data[model][lo:hi], Rand[lo:hi], Years[lo:hi], np.array([1880.0, 1892.0, 1931.0, 1970.0,  1983.0, 1995.0,2012.0]))
    print "Linear Rand",a.evaluate(data[model][lo:hi], Rand[lo:hi], Years[lo:hi], np.array([1880.0, 1892.0, 1931.0, 1968,  1983.0, 1995.0, 2012.0]))
    print "Linear Rand",a.evaluate(data[model][lo:hi], Rand[lo:hi], Years[lo:hi], np.array([1880.0, 1892.0, 1931.0, 1970.0,  1995.0, 2012.0]))
 
    #raise Exception

    a=AIC(data[model][lo:hi], Years[lo:hi], Years[lo:hi], np.array([1880.0, 1892.0, 1931.0, 1968.0, 1995.0, 2012.0]))
    print "Linear Years", a.value()
    print "Linear Years", a.evaluate(data[model][lo:hi], Years[lo:hi], Years[lo:hi], np.array([1880.0, 1892.0, 1931.0, 1970.0,  1983.0, 1995.0,2012.0]))
    print "Linear Years", a.evaluate(data[model][lo:hi], Years[lo:hi], Years[lo:hi], np.array([1880.0, 1892.0, 1931.0, 1968,  1983.0, 1995.0, 2012.0]))
    print "Linear Years", a.evaluate(data[model][lo:hi], Years[lo:hi], Years[lo:hi], np.array([1880.0, 1892.0, 1931.0, 1970.0,  1995.0,2012.0]))

#[1880.0, 1892.0, 1931.0, 1968.0, 1970.0, 1995.0, 2012.0]
