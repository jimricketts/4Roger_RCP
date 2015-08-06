#!/usr/bin/env python
#support for the 
import numpy as np
import scipy
import os
import sys
from scipy import stats
from scipy.stats import t
import math
import string
import glob
import subprocess as sp
#import Nio
#import Ngl
#import gcm
import globalw
import numpy.ma as ma
import copy
import regress
import shutil
import matplotlib 
#matplotlib.use("Agg") # or whichever backend you wish to use ..
import matplotlib.pyplot as plt
from math import log
import random
#import pbd #<<<<<<<<<<<<<<<debugger
SVNRevision="$Revision: 307 $"
import linecache
#import sys

def ExceptionInfo():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    return 'EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj)
   
    
class BivariateException(Exception):
  pass
  
#class BivariateException(Exception):
#    def __init__(self, value):
#        self.value = value
#    def __str__(self):
#        return repr(self.value)
  def __str__(self):
    return "BivariateException: "+str(self.message)

def Pr2(Ti1, N1): #Alternative model obtained 26 august 
  Ti = max(2.0, min(14.0, Ti1))
  N = max (5, min(N1, 100))
  return max(0.000001,0.995617028964232*1.00034593706459**N*0.624183483006978**Ti*log(15.2050131939315 + 0.60200250374432*N - 2.4951612324163*Ti))

def Pr(Ti, N):
  '''
  Significance level (Jones and Vives 2005)
n (years) 0.25 0.10 0.05 0.01
10        4.7  6.0   6.8  7.9
15        4.9  6.5   7.4  9.3
20        5.0  6.7   7.8  9.8
30        5.3  7.0   8.2  10.7
40        5.4  7.3   8.7  11.6
70        5.9  7.9   9.3  12.2
100       6.0  7.9   9.3  12.5

Eureqa Version:  0.99.8 Beta (build 4023)
Revision: 
f39f6ebb4842733696cea08cd128db7a95beed8f
Architecture: 
64 bit
Operating System: 
Linux 
Language: 
English
Country: 
UnitedStates
CPU Name: 
Intel(R) Xeon(R) CPU 5140 @ 2.33GHz
CPU Cores: 
4
UUID: 
2486bba2-4e34-4741-beee-58c1b393e39f
Access ID: 


Pr = 0.0469591507344088*Ti + 0.00530855792450758*N + 4.27735087880932/Ti + -28.6534878246636/N^3 + 1.56788437630966e-7*N^2*Ti^2 - 0.883614188389141 - 0.000558194019548398*N*Ti - 1.13551662914624e-7*N^3

Alternative model obtained 26 august 
Pr = 0.995617028964232*1.00034593706459^N*0.624183483006978^Ti*log(15.2050131939315 + 0.60200250374432*N - 2.4951612324163*Ti)

  '''
  return max(0.000001,0.0469591507344088*Ti + 0.00530855792450758*N + 4.27735087880932/Ti + -28.6534878246636/N**3 + 1.56788437630966e-7*N**2*Ti**2 - 0.883614188389141 - 0.000558194019548398*N*Ti - 1.13551662914624e-7*N**3)

def critTi(Pr, N1, lowerlim=10, upperlim=200):

  #print "critTi", Pr, type(Pr), N, type(N)
  if np.shape(N1) == ():
    N = max(lowerlim, N1)
    N = min(N1, upperlim)
    p1 = N -0.0452658289957312/(N*Pr - 8.48780080143555*Pr)
    if p1 <= 0.0 :
      return 4.29938163981431 - 3.68238858146661*Pr
    else:
    #"Ti = 4.29938163981431 - 3.68238858146661*Pr - 0.405718682094657*log(Pr)*log(N + -0.0452658289957312/(N*Pr - 8.48780080143555*Pr))"
      return 4.29938163981431 - 3.68238858146661*Pr - 0.405718682094657*log(Pr)*log(p1)
  else:
    N = np.where(N1 < lowerlim, lowerlim, np.where( N1 > upperlim, upperlim, N1))
    p1 = N -0.0452658289957312/(N*Pr - 8.48780080143555*Pr)
    t2=np.where(p1 <= 0.0, 4.29938163981431 - 3.68238858146661*Pr, 4.29938163981431 - 3.68238858146661*Pr - 0.405718682094657*np.log(Pr)*np.log(p1))
    return np.min(8.0, t2)

def PrEst(Ti, N):
  if N <= 100: 
    return Pr(Ti, N)
  else: #out of range. Use a bisection method
    crit=critTi(0.000001, N)
    if Ti > crit:
      #print crit
      return 0.000001
    else:
      #bisection
      low=0.000001
      hi = 1.0 - 0.000001
      mid = (hi + low)/2.0
      midTi = critTi(mid, N)
      loop = True
      while loop:
       # print hi, low, mid, midTi, hi > low
        if midTi < Ti:
          hi = mid
        elif midTi > Ti:
          low = mid
        else:
          hi = low
        mid = (hi + low)/2.0
        midTi = critTi(mid, N)
        loop = abs(low - hi) > 0.000001
    return mid

def _chk_asarray(a, axis):
    if axis is None:
        a = np.ravel(a)
        outaxis = 0
    else:
        a = np.asarray(a)
        outaxis = axis
    return a, outaxis


def nanstd(x, msg=None, axis=0, bias=False):
    """
Compute the standard deviation over the given axis, ignoring nans.

Parameters
----------
x : array_like
Input array.
axis : int or None, optional
Axis along which the standard deviation is computed. Default is 0.
If None, compute over the whole array `x`.
bias : bool, optional
If True, the biased (normalized by N) definition is used. If False
(default), the unbiased definition is used.

Returns
-------
s : float
The standard deviation.

See Also
--------
nanmean, nanmedian

Examples
--------
>>> from scipy import stats
>>> a = np.arange(10, dtype=float)
>>> a[1:3] = np.nan
>>> np.std(a)
nan
>>> stats.nanstd(a)
2.9154759474226504
>>> stats.nanstd(a.reshape(2, 5), axis=1)
array([ 2.0817, 1.5811])
>>> stats.nanstd(a.reshape(2, 5), axis=None)
2.9154759474226504

"""
    x, axis = _chk_asarray(x, axis)
    x = np.array(x.copy(), dtype = np.float64) #<cast to 64 bit reduces a problem due to rounding. But algorithms is naive.
    Norig = x.shape[axis]

    mask = np.isnan(x)
    Nnan = np.sum(mask, axis) * 1.0
    n = Norig - Nnan

    x[mask] = 0.0
    m1 = np.sum(x, axis) / n
#    pbd.set_trace()
    if axis:
        d = x - np.expand_dims(m1, axis)
        #if msg != None: print "DATASTD is axis", msg, n, m1, 
    else:
        d = x - m1
        #if msg != None: print "DATASTD not axis",msg, n, m1,
    d *= d
    
    m2 = np.sum(d, axis) - m1 * m1 * Nnan #<<<<<<<<<<<<<<<<<<<<<<<<<This is what breakes
    #if msg != None: print Nnan, m2, bias,
    if bias:
        m2c = m2 / n
    else:
        m2c = m2 / (n - 1.0)
   # if msg != None: print m2c, np.sqrt(m2c)
    return np.sqrt(m2c)


def resample_break(testdata, N=30, withmode=True):
  brks=[]
  tis=[]
  shifts=[]
  try:
    for i in range(N):
      step=0
      controldata = np.array([random.random() for y in testdata])
      step=1
      bv=bivariate(testdata,controldata, anomalise=False, pr=0.01)
      step=2
      tis.append(bv.maxTi())
      step=3
      brks.append(bv.maxIndexTi())
      step=4
      shifts.append(bv.stepChange())
      step=5
  except Exception as e:
    
    print str(e), step
    raise BivariateException(str(e)+"\n"+ExceptionInfo())
  try:
    if withmode:
      yearfreqs=np.bincount(brks)
      #print brks,len(yearfreqs),yearfreqs,range(int(datayears[0])-1,int(datayears[-1])-1)
      first=second=firstval=secondval=0
      for i in range(int(min(brks))-1,int(max(brks))+1):
       # print i
        if yearfreqs[i] >firstval:
          second=first
          first=i
          firstval=yearfreqs[first]
          secondval=yearfreqs[second]
        elif yearfreqs[i] >secondval:
          second=i
          secondval=yearfreqs[second]
      #now we also will need the mean of the Tis and shifts for each of the modal values
      brks=np.array(brks)
      mask1=np.where(brks==first) 
      timean1=np.mean(np.array(tis)[mask1])
      shmean1=np.mean(np.array(shifts)[mask1])
      if second == 0:
        return stats.norm.fit(brks[mask1]), stats.norm.fit(np.array(tis)[mask1]),stats.norm.fit(np.array(shifts)[mask1]), [(first, float(firstval)/N), None, None,(timean1, shmean1), None]
      else:
        mask2=np.where(brks==second) 
        timean2=np.mean(np.array(tis)[mask2])
        shmean2=np.mean(np.array(shifts)[mask2])
        return stats.norm.fit(brks[mask1]), stats.norm.fit(np.array(tis)[mask1]),stats.norm.fit(np.array(shifts)[mask1]), [(first, float(firstval)/float(N)), (second, float(secondval)/float(N)), yearfreqs[-len(testdata):],(timean1, shmean1), (timean2, shmean2)]
  #    if second == 0:
  #      return stats.norm.fit(brks), stats.norm.fit(tis),stats.norm.fit(shifts), [(first, float(firstval)/N), None, None,(timean1, shmean1), None]
  #    else:
  #      mask2=np.where(brks==second) 
  #      timean2=np.mean(np.array(tis)[mask2])
  #      shmean2=np.mean(np.array(shifts)[mask2])
  #      return stats.norm.fit(brks), stats.norm.fit(tis),stats.norm.fit(shifts), [(first, float(firstval)/float(N)), (second, float(secondval)/float(N)), yearfreqs[-len(datayears):],(timean1, shmean1), (timean2, shmean2)]
    else:
      return stats.norm.fit(brks), stats.norm.fit(tis),stats.norm.fit(shifts)  
  except Exception as e:
    raise BivariateException(str(e)+"\n"+ExceptionInfo())


        
dumpme = False
class bivariate(object):
    class breakclassification(object):
      def __init__(self, bv):
        self.__str="str()" 
        
    def __init__(self, rawdata, xs, anomalise=True, step=1, averagesteps=False, critical=None, pr= None, window = 5, constantsxy=True):
        #print "entering init"
        #save parameters
        self.anomalise = anomalise
        self.step=step, 
        self.averagesteps=averagesteps
        self.criticalin=critical
        self.pr= pr
        self.window =window
        self.constantsxy=constantsxy
        self.__storedStep=None
    #0. if monthly data is used removed the mean monthly values (and the annual cycle)
        
        if anomalise:
            shape=np.shape(rawdata)
            shape=(12, shape[1], shape[2])
            monthly=np.zeros(shape)
            data=np.copy(np.array(rawdata))
            for i in range(12):
                monthly[i] = np.mean(data[i::12], axis=0)
            for i in range(len(data)):
                data[i] = data[i] - monthly[ i % 12]
        else:
            data=np.copy(rawdata)
    #0.1 if stepsite of anything greater than 1 is specified deal with eaxtraction / averaging across that step
        if critical == None:
          if pr == None:
            self.__critical = critTi(0.01, len(xs))
          else:
            self.__critical = critTi(pr, len(xs))
        else:
          self.__critical = critical
        self.__window = window
        if step <1:
            print "bivariate ignoring step size (must be >= 1)"
        elif step > 1:
            if averagesteps:
                data = np.array([np.mean(data[i*step:(i+1)*step],axis=0) for i in range(int(len(data) / step))])
                xs=np.array([np.mean(xs[i*step:(i+1)*step],axis=0) for i in range(int(len(data) / step))])
                
            else:
                data = np.array(data)[::step]
                xs = np.array(xs)[::step]
    #0.2 save the unnormalised data raw data for later reuse
        self.__rawadata=np.copy(data)
        self.__rawxs = np.copy(xs)
    
    #1. normalisation of data
        #print "normalise"
        datastd=nanstd(data,"DATA",axis=0,bias=False)
        datamean=stats.nanmean(data,axis=0)
        xmean=stats.nanmean(xs,axis=0)
        xstd=nanstd(xs, "XS",axis=0, bias=False)
        self.__xi=(xs-xmean)/xstd
        self.__yi=(data-datamean)/datastd
        if dumpme:
          for i in range(len(data)):
            print "DUMP_1", i, data[i]

    #2. running averages
        #copy the yis to a another array, similarly xs
        self.__Yi=np.copy(self.__yi)
        self.__Xi=np.copy(self.__xi)
        #now create running sums of the arrays
        for i in range(1,len(self.__Xi)):
            self.__Yi[i] +=self.__Yi[i-1]
            self.__Xi[i] +=self.__Xi[i-1]
        #turn these into averages        
        for i in range(1,len(self.__Xi)):
            self.__Yi[i] /=(i+1)
            self.__Xi[i] /=(i+1)
            
    #3. Create Sxy  1-jun-2014 modified to test my theory that Sxy is not intended to vary
        self.__constantsxy=constantsxy
        if constantsxy:
          sxy=np.sum(self.__yi * self.__xi, axis=0)
          self.__Sxy=np.array([sxy for i in range(len(self.__yi))])
        else:
          self.__Sxy=np.copy(self.__yi)
          for i in range(len(self.__xi)):
              self.__Sxy[i]=self.__Sxy[i]*self.__xi[i]
          for i in range(1,len(self.__xi)):
              self.__Sxy[i] += self.__Sxy[i-1]  
          
            
    #4. some other useful things
        self.__N = len(self.__xi)       
        ns=np.array([np.ones(np.shape(self.__xi[0])) * t for t in range(len(self.__xi))]) + 1
        self.__ns = ns
        self.__datastd = datastd
        self.__datamean = datamean
        self.__xstd = xstd
        self.__xmean = xmean
        if dumpme: print "DATASTD_1", datastd, datamean, self.__N, ns, self.__yi[:]
        
    #5. Create Fi
        try:
          self.__Fi=np.copy(self.__N - (self.__Xi * self.__Xi * ns * self.__N )/(self.__N - ns))
        except:
          pass
            
    #6. Create Di
        self.__Di = np.copy(self.__yi)
        for i in range(len(self.__xi)):
          try:
            self.__Di[i]=((self.__Sxy[i]*self.__Xi[i])-(self.__N*self.__Yi[i]))*self.__N/((self.__N-ns[i])*self.__Fi[i])
          except:
            pass
    #7  Create Ti
        self.__Ti = np.copy(self.__yi)
        for i in range(len(self.__xi)):
          try:
            self.__Ti[i] = ns[i] * (self.__N - ns[i]) * self.__Di[i] * self.__Di[i] * self.__Fi[i]/(self.__N * self.__N - self.__Sxy[i] * self.__Sxy[i])
          except:
            pass
        #self.__stats=regress.analysed_regress(data,xs)
   
    #8 Locate the maxTi 
   #JHR 28/5/2014 deal will all nuls case
        if not np.isnan(self.__Ti[1:-1]).all():
          self.__maxTi= np.nanmax(self.__Ti[1:-1],axis=0)    
          #self.__maxIndexTi = np.nanargmax(self.__Ti[1:-1], axis=0)
          #JHR 20aug2014 This is incorrect because it is returning Toi indexed once early
          self.__maxIndexTi = np.nanargmax(self.__Ti[1:-1], axis=0) + 1
        else:
          self.__maxTi= np.NaN
          self.__maxIndexTi = 0
         
        
        
    def stats(self):
        return self.__stats
      
    def N(self):
        return self.__N
      
    def ns(self):
        return self.__ns
    
    def Xi(self):
        return self.__Xi
    
    def Yi(self):
        return self.__Yi

    def xi(self):
        return self.__xi
      
    def yi(self):
        return self.__yi
      
    def Sxy(self):
        return self.__Sxy       
        
    def Fi(self):
        return self.__Fi[:-1]  

    def Di(self):
        return self.__Di[:-1]  

    def Ti(self):
        return self.__Ti[:-1]       

    def maxTi(self, Array=False):
        if Array and np.shape(self.__maxTi) == ():
          return np.array([self.__maxTi])
        else: 
          return self.__maxTi     
        
    def stepChange(self):
        if self.__storedStep == None:
          return np.max(self.__Di[self.__maxIndexTi]) * self.__datastd
        else:
          return self.__storedStep
        
    def maxIndexTi(self):
        return self.__maxIndexTi  
      
    def datastd(self):
        return self.__datastd
      
    def datamean(self):
        return self.__datamean
      
    def xstd(self):
        return self.__xstd
      
    def xmean(self):
        return self.__xmean
      
    def rawxs(self):
      return self.__rawxs
    
    def rawadata(self):
      return self.__rawadata
    
    #JHR 22/5/2014 the chnage in mean at nominated point or at -I don't think this is correct
    def delta(self, PtNumber=None):
      if PtNumber==None:
        ix = self.__maxIndexTi
      else:
        ix = PtNumber
      return self.__Di[ix] *self.__Yi[ix] 
    
    def reinit(self, low, high, pr, bymodes, debug=False): #assume 2D array
        if np.shape(low) != np.shape(high) or np.shape(low) != np.shape(self.maxTi()):
            raise BivariateException ("reinit() index shape clash, shape(__Yi[0])=%s shape(low)=%s shape(high)=%s MaxTi()=%s low=%s high=%s" % 
                  (str(np.shape(self.__Yi[0])), str(np.shape(low)), str(np.shape(high)), str(self.maxTi()), str(low), str(high)))
 
            sys.exit()
            #refurmulate here (data already anomalised)
        try:
          tempbv=bivariate(self.rawadata()[int(low):int(high)], self.rawxs()[int(low):int(high)], anomalise=False, step=1, averagesteps=self.averagesteps,critical=self.criticalin, pr=self.pr,window=self.window,constantsxy=self.constantsxy)
          if bymodes:
            brks_, tis_, shifts_, modes_ = resample_break(self.rawadata()[int(low):int(high)], N=100, withmode=True)
            try:
              tempbv._bivariate__maxIndexTi=int(brks_[0])
              tempbv._bivariate__maxTi=tis_[0]
              tempbv._bivariate__storedStep=shifts_[0]
            except:
              tempbv._bivariate__maxIndexTi=None
              tempbv._bivariate__maxTi=0.0
              tempbv._bivariate__storedStep=0.0
              pass
        except Exception as e:
          einfo=ExceptionInfo()
          print einfo
          #raise BivariateException(str(e)+"\n"+einfo)
          tempbv = None
          pass          
        return tempbv
    '''
    def maxnTi(self):
        return self.__tempbv.maxTi()     
        
    def maxIndexnTi(self):
        return self.__tempbv.maxIndexTi()
      
    def nSxy(self):
        return self.__tempbv.Sxy()
      
    def nTi(self):
        return self.__tempbv.Ti()
    
    def nFi(self):
        return self.__tempbv.Fi()
      
    def nDi(self):
        return self.__tempbv.Di()
      
    def nNs(self):
        return self.__tempbv.Ns() 
      
    def nns(self):
        return self.__tempbv.ns()
         
    def mask(self):
        return self.__tempbv.mask()
      
    def summask(self):
        return self.__tempbv.__summask
         
    def nXi(self):
        return self.__tempbv.Xi()
      
    def nYi(self):
        return self.__tempbv.Yi()
      
    def ndatastd(self):
        return self.__tempbv.datastd()
      
    def ndatamean(self):
        return self.__tempbv.datamean()
      
    def nxstd(self):
        return self.__tempbv.xstd()
      
    def nxmean(self):
        return self.__tempbv.xmean()

    def nxi(self):
        return self.__tempbv.xi()
      
    def nyi(self):
        return self.__tempbv.yi()
    '''  
    def allPoints(self, pr, debug=False, allti=False, withshifts=False, bymodes=True):
      #already initialised so we use the maxIndexTi
      #and just return lists of candidate points and their t values
      #get the mid points
      if bymodes:
        #then we discard all to date and start again.
        brks_, tis_, shifts_, modes_ = resample_break(self.rawadata()[:], N=100, withmode=True)
        ap_MaxTis = [np.where(tis_[0] >= self.__critical, tis_[0], np.NaN)]
        ap_MaxIndexes = [np.where(np.isnan(ap_MaxTis[0]), np.NaN, brks_[0])] #had been failing to dereference MaxTis JHR 5/3/2014
        ap_lows=[np.zeros(np.shape(self.__maxTi))]
        ap_highs=[np.array(np.ones(np.shape(self.__maxTi))*len(self.xi()[:]))]
        if withshifts: 
          shft = shifts_[0]
          ap_shifts=[np.where(np.isnan(ap_MaxTis[0]), np.NaN, shft)]
      
      else:
        ap_MaxTis = [np.where(self.__maxTi >= self.__critical, self.__maxTi, np.NaN)]
        #print MaxTis, self.__maxTi, self.__Ti
        #sys.exit()
        ap_MaxIndexes = [np.where(np.isnan(ap_MaxTis[0]), np.NaN, self.__maxIndexTi)] #had been failing to dereference MaxTis JHR 5/3/2014
        ap_lows=[np.zeros(np.shape(self.__maxTi))]
        ap_highs=[np.array(np.ones(np.shape(self.__maxTi))*len(self.xi()[:]))]
        if withshifts: 
          shft = self.__Di[self.__maxIndexTi] * self.__datastd
          ap_shifts=[np.where(np.isnan(ap_MaxTis[0]), np.NaN, shft)]
#>>>The change    self.__nDi[self.__maxIndexnTi]) * self.__ndatastd   
      
      #now define an embedded rucursive function that 'cheats' using globally scoped variables. I may regret this. 29/2/2014 
      
      if np.isnan(ap_MaxTis).all():
        ap_MaxTis = []
        ap_MaxIndexes=[]
        ap_lows=[] 
        ap_highs=[]
        if withshifts:
          return [],[],([],[]),[]
        else:
          return [],[],([],[])
      
      #capture the new Ti values - set up the array here and update on exit of each recursion JHR 9/6/2014
      self.__recursiveTi=None
      def recurse(lows, highs, mids, pr, trace=False):
        #define a recursive search to search between low and mid, and between mid and high
        #skip if there is nothing interesting
        try:
          if trace: print "Enter Recursion " , lows, highs, mids
          if debug: print lows, highs, "\n",mids
          if lows != None:
            if debug: print "LOW ",
            try:
             # print np.shape(lows), np.shape(highs)
              try:
   ##fail 1             
                lbv=self.reinit(np.array(lows), np.array(mids), pr, bymodes)
              except Exception as e:
                print bivariate.ExceptionInfo()
                print "RIHGHT", np.array(lows), np.array(mids), pr, ExceptionInfo()
                raise "RIHGHT " + str(np.array(lows)) + " " + str(np.array(mids)) + " " + str(e)
              if dumpme: print "MEANNSES LOW", lows, mids, lbv.datamean(), lbv.xmean(), lbv.datastd(), lbv.xstd()
              if debug:
                try:
                  print "reinit nTi", lows, highs, mids, lbv.Ti(), lbv.ns()
                  for k in range(lows,highs):
                    print k, lbv.ns()[k], lbv.Sxy()[k], lbv.Xi()[k],lbv.Fi()[k], lbv.Di()[k], lbv.Ti()[k], lbv.xi()[k], lbv.yi()[k], lbv.N(), lbv.Ns()
                except:
                  pass
            except Exception as e:
              print ExceptionInfo()
              raise BivariateException ("Low error %s %s " % (str((lows,highs)), str(mids)) + str(e))
            lMaxTis=lbv.maxTi()
            if debug: print lMaxTis
                  
            #replaed this lowMaxTis=np.where(lMaxTis >= self.__critical and mids > lows and mids < highs and highs - lows >= self.__window, lMaxTis, np.NaN)
            if debug: print type(lbv.maxIndexTi())
            if lbv.maxIndexTi() == None:
              bvmi = lows
            else:
              bvmi = np.where(np.isnan(lbv.maxIndexTi()), lows, lows+lbv.maxIndexTi()) ##JHR
            lowMaxTis=np.where(lMaxTis < lbv.critical(), np.NaN, 
                               np.where(mids <= lows, np.NaN,
                                        np.where(mids >= highs, np.NaN,
                                                np.where(abs(bvmi-lows) < self.__window, np.NaN, lMaxTis))))
            
            lowMaxIndexnTis=np.where(np.isnan(lowMaxTis), np.NaN, bvmi)
            if debug: print "lowMaxIndexnTis",lowMaxIndexnTis
            ##JHR 7/5/14 this is the place to save a potential list of break points.
            ##JHR 26/8/14 also save shift values
            if not np.isnan(lowMaxTis).all():
              ap_lows.append(np.copy(lows))
              ap_highs.append(np.copy(mids))
              ap_MaxTis.append(np.copy(lowMaxTis))
              ap_MaxIndexes.append(np.copy(lowMaxIndexnTis))
              if withshifts: ap_shifts.append(np.where(np.isnan(lowMaxTis), np.NaN, lbv.Di()[lbv.maxIndexTi()]) * lbv.datastd())
              #capture the new Ti values
              #JHR No real use for this at present, just append to list
              if self.__recursiveTi == None:
                self.__recursiveTi = []
              self.__recursiveTi.append(lbv.Ti())
#                self.__recursiveTi=np.zeros(np.shape(lbv.Ti()))
#              self.__recursiveTi=np.where(np.isnan(lbv.Ti()), self.__recursiveTi, lbv.Ti())
            else:
              lowMaxTis = None
              lowMaxIndexnTis = None
            
          if highs != None:
            if debug: print "HIGH",
            try:
              hbv=self.reinit(np.array(mids), np.array(highs), pr,bymodes )
              #bv=self.reinit(np.array(mids)+1, np.array(highs), pr ) #I think this is correct but it crashed JHR 20auf14
              if debug: 
                try:
                  print "reinit nTi", lows, highs, mids, hbv.Ti(), hbv.ns()
                  for k in range(lows,highs):
                    print k, hbv.ns()[k], hbv.Sxy()[k], hbv.Xi()[k],hbv.Fi()[k], hbv.Di()[k], hbv.Ti()[k]
                except:
                  pass
            except BivariateException as e:
              print ExceptionInfo()
              raise BivariateException ("High error %s %s " % (str((lows,highs)), str(mids)) + str(e))
            hMaxTis=hbv.maxTi()
            if hbv.maxIndexTi() == None:
              bvmi = highs
            else:
              bvmi = np.where(np.isnan(hbv.maxIndexTi()), highs,mids + hbv.maxIndexTi())
            if debug: print hMaxTis
            #replace highMaxTis=np.where(hMaxTis >= self.__critical and mids > lows and mids < highs and highs - lows >= self.__window, hMaxTis, np.NaN)
            #22/05/2014 JHR the line below used lMaxTis instead of hMaxTis - problems introduced as commented above.
            #7/3/2015 JHR this isalso incorrect.
#            print "hMaxTis < hbv.critical()",hMaxTis < hbv.critical()
#            print "abs(bvmi-highs) < self.__window", abs(bvmi-highs) < self.__window
#            print " np.where(abs(bvmi-mids) >=self.__window", abs(bvmi-mids) >=self.__window
            highMaxTis=np.where(hMaxTis < hbv.critical(), np.NaN, 
                               np.where(mids <= lows, np.NaN,
                                        np.where(mids >= highs, np.NaN,
                                                np.where(abs(bvmi-highs) >= self.__window, 
                                                         np.where(abs(bvmi-mids) >=self.__window, hMaxTis, np.NaN), np.NaN))))
  
            highMaxIndexnTis=np.where(np.isnan(highMaxTis) , np.NaN, bvmi) #JHR 9/3/2015 this was incorrect ater mods
            ##JHR 7/5/14 this is the place to save a potential list of break points.
            if not np.isnan(highMaxTis).all():
              ap_lows.append(np.copy(mids))
              ap_highs.append(np.copy(highs))
              ap_MaxTis.append(np.copy(highMaxTis))
              ap_MaxIndexes.append(np.copy(highMaxIndexnTis))
              if withshifts: ap_shifts.append(np.where(np.isnan(highMaxTis), np.NaN, hbv.stepChange()))
              #capture the new Ti values
              #JHR No real use for this at present, just append to list
              if self.__recursiveTi == None:
                self.__recursiveTi = []
              self.__recursiveTi.append(hbv.Ti())
#                self.__recursiveTi=np.zeros(np.shape(hbv.Ti()))
#              
#              #print np.shape(bv.nTi()), np.shape(bv.__recursiveTi)
#              self.__recursiveTi=np.where(np.isnan(hbv.Ti()), self.__recursiveTi, hbv.Ti())
            else:
              highMaxTis = None
              highMaxIndexnTis = None
          #do we go on?
  #==============================================================================
  #         if np.isnan(lowMaxTis).all(): 
  #           lowMaxTis = None
  #           lowMaxIndexnTis = None
  #         if np.isnan(highMaxTis).all():
  #           highMaxTis = None
  #           highMaxIndexnTis = None
  #==============================================================================
          #save results if needed - JHR 7/5/14 this is wwrong and I'm sure I fixed it. Not here apparently
  #==============================================================================
  #         if lowMaxTis != None or highMaxTis != None:
  #           ap_lows.append(np.copy(lowMaxIndexnTis))
  #           ap_highs.append(np.copy(highMaxIndexnTis))
  #           ap_MaxTis.append(np.copy(lowMaxTis))
  #           ap_MaxIndexes.append(np.copy(lowMaxIndexnTis))
  #==============================================================================
            #now recursivly explore the subsections
            #clean up the bivariates first
            lbv = None
            hbv = None
            if lowMaxTis != None:
              try:
                if ((mids-lows) != 0.0).any():
                  recurse(np.copy(lows), np.copy(mids), lowMaxIndexnTis, pr)
              except Exception as e:
                print ExceptionInfo()
                #recurse(np.copy(lows), np.copy(mids), lowMaxIndexnTis, pr)
                print str(e), " low recurse failed", np.copy(lows), np.copy(mids), lowMaxIndexnTis, pr
                pass
            if highMaxTis != None:
              try:
                if ((highs - mids) != 0).any():
                  recurse(np.copy(mids), np.copy(highs), highMaxIndexnTis, pr)
              except Exception as e:
                print ExceptionInfo()
                #recurse(np.copy(mids), np.copy(highs), highMaxIndexnTis, pr)
                print str(e), " high recurse failed", np.copy(mids), np.copy(highs), highMaxIndexnTis, pr
                pass
          #else:
            #ap_lows.append(np.array(None))
            #ap_highs.append(np.array(None))
          if trace: print "Exit Recursion ", lows, "-> ", mids, lowMaxTis, lowMaxIndexnTis, mids, "->", highs, highMaxTis, highMaxIndexnTis
          #recursions += 1
          #and finally just return
          #return recursions
        except Exception as e:
          print ExceptionInfo()
          raise e
          
        return None  ##end def recurse
        
      #if debug: 
      #print lows, highs, MaxIndexes
      recurse(ap_lows[0], ap_highs[0], ap_MaxIndexes[0], pr)
      if len(ap_MaxTis) != len(ap_MaxIndexes):
        raise BivariateException("len(ap_MaxTis) != len(ap_MaxIndexes) (%d != %d)" % (len(ap_MaxTis), len(ap_MaxIndexes)))
        
      #capture the new Ti values here
        
      if allti:
        return ap_MaxTis, ap_MaxIndexes, (ap_lows, ap_highs), self.__recursiveTi
      elif withshifts:
        return ap_MaxTis, ap_MaxIndexes, (ap_lows, ap_highs), ap_shifts      
      else:
        return ap_MaxTis, ap_MaxIndexes, (ap_lows, ap_highs)
    
    def window(self):
      return self.__window
    
    def critical(self):
      return self.__critical
    
    def ncritical(self):
      return self.__ncritical
      
    def allTis(self):
      if self.__recursiveTi == None:
        return self.Ti()
      return self.__recursiveTi[-1]
      
class ti_accumulator(object):
  def __init__(self):
    self.__tilist=[]
  
  def append(self, bv):
    self.__tilist.append(bv.allTis())
    
  def classify(self, ap_MaxIndex, ap_low=None, ap_high=None):
    print "TILIST",self.__tilist
    dist=(np.array([range(len(t)) for t in self.__tilist],dtype=np.float32) -ap_MaxIndex).flatten() 
    y=np.array(self.__tilist,dtype=np.float32).flatten()
    
    return dist, y
    

#def breakclass(bv, ap_MaxTis, ap_MaxIndexes, ap_lows, ap_highs): 
#  '''
#  classify the break around Ti0 and return a breakclassification object
#      
#  '''
    
  
if __name__ == "__main__":
    import sys
    import somedata
    import victoria
    import random
    Years = np.array([victoria.VictorianAnnual['data'][i][0] for i in range(len(victoria.VictorianAnnual['data']))])
    Tmean = np.array([victoria.VictorianAnnual['data'][i][1] for i in range(len(victoria.VictorianAnnual['data']))])
    RainAnom = np.array([victoria.VictorianAnnual['data'][i][2] for i in range(len(victoria.VictorianAnnual['data']))])
    Tmax = np.array([victoria.VictorianAnnual['data'][i][3] for i in range(len(victoria.VictorianAnnual['data']))])
    Tmin = np.array([victoria.VictorianAnnual['data'][i][4] for i in range(len(victoria.VictorianAnnual['data']))])
    Diurnal = np.array([victoria.VictorianAnnual['data'][i][5] for i in range(len(victoria.VictorianAnnual['data']))])
    xs = np.array([random.random() for i in range(len(victoria.VictorianAnnual['data']))])
    #print Years, Diurnal
    bv=bivariate(Tmean,Diurnal, pr = 0.05, anomalise=False)
    ap = bv.allPoints(0.05,debug=False)
    print "AP",ap
#    sys.exit()
#    
    tiac=ti_accumulator()
    for i in range(10):
      xs = np.array([random.random() for i in range(len(victoria.VictorianAnnual['data']))])
      bv=bivariate(Tmean[:88],xs[:88], pr = 0.05, anomalise=False)
      ap=bv.allPoints(pr = 0.05)
      print "TAP",ap
      tiac.append(bv)
    x, y =tiac.classify(59)
    for i in range(len(x)):
      print x[i],y[i]
    sys.exit()    
    #MyData=np.array(somedata.temps)[:]
    #column = 1
    #ys=np.array([d[column] for d in MyData], np.float32)    
    #xs=np.array([np.random.randn() for i in range(len(MyData))], np.float32)    
    ##xs=np.array([random.random()for d in MyData], np.float32)    
    ##print somedata.Header[i], somedata.Temps[i]
    #i = len(MyData)             

    #bv=bivariate(ys[:], xs[:], anomalise=False)        
    #Ti=bv.maxTi()
    #i=bv.maxIndexTi()
    #print "Full",i, Ti, somedata.Header[column], MyData[i], "\n",bv.Sxy()
    
    #bv=bivariate(ys[:87], xs[:87], anomalise=False)        
    #Ti=bv.maxTi()
    #i=bv.maxIndexTi()
    #print "BP 1",i, Ti, somedata.Header[column], MyData[i], "\n",bv.Sxy()

    #bv=bivariate(ys[:], xs[:], anomalise=False)        
    #bvv=bv.reinit(0,87)        
    #Ti=bv.maxnTi()
    #i=bv.maxIndexnTi()
    #print "New",i, Ti, somedata.Header[column], MyData[i], "\n",bvv.nSxy()
    #sys.exit()
#AN INTERPRETER   
    xData=[0.08,0.03,0.05,0.02,0.00,0.10,0.05,0.06,0.03,0.02,0.06,0.08,0.01,0.06,0.09,0.05,0.04,0.05,0.05,0.04,0.08,0.05,0.03]

    Ydata=np.zeros((23,2,2))
    j=0
    for data in [1.64,1.64,1.61,1.82,1.84,1.76,1.66,1.67,1.81,1.91,1.92,1.92,1.83,1.89,2.04,1.98,2.13,2.06,1.99,2.08,2.11,2.16,2.19]:
       Ydata[j][0][0]=data 
       Ydata[j][1][0]=data 
       Ydata[j][0][1]=data 
       Ydata[j][1][1]=data 
       j +=1
#    print bivariate(Ydata, xData, anomalise=False).maxIndexTi()
    import somedata
    import random
    DoRand=False
    
    MyData=np.array(somedata.temps)[:]
    try:
      with open('MyXdata.txt','r') as xf:
        MyxData=np.array([float(line) for line in xf.readlines()])
        #print MyxData
    except:
        MyxData=np.array([random.random()for d in MyData], np.float32)   
        try:
          with open('MyXdata.txt','w') as xf:
            for dat in MyxData:
              xf.write("%s\n" % (str(dat),))
        except Exception as e:
          raise BivariateException('error with Myxdata.txt %s ' %(e))
    #column=1
    #print "Now"
    #for column in [1]:
        #ys=np.array([d[column] for d in MyData], np.float32)    
        #xs=np.array([np.random.randn() for i in range(len(MyData))], np.float32)    
        ##xs=np.array([random.random()for d in MyData], np.float32)    
        ##print somedata.Header[i], somedata.Temps[i]
        #i = len(MyData)
        #bv=bivariate(ys[:], xs[:], anomalise=False)
        ##print "xi yi\n",bv.xi(), "\n",bv.yi()
##        bvv=bv.reinit(np.zeros(np.shape(bv.maxIndexTi())), bv.maxIndexTi())       
        #bvv=bv.reinit(np.zeros(np.shape(bv.maxIndexTi())), len(xs[:])* np.ones(np.shape(bv.maxIndexTi())),debug=False)       
        #Ti=bv.maxnTi()
        #i=bv.maxIndexnTi()
        #print i, Ti, somedata.Header[column], MyData[i]
        #print "Tis, nTi then Ti"
        #print bvv.nTi(), "\n",bv.Ti(),"\n",bvv.nTi()[:122]/bv.Ti()[:122]
        ##print "xi yi 2\n",bv.xi(), "\n",bv.yi()
##        print bv.xi(), bv.yi()
    #sys.exit()        
#######################################################################################################################   

#testing the recursive stuff
#multiples
    tasD={}
    if len(sys.argv) <=1:
      sys.argv.append("C:\Users\\s4493222\\Documents\\abrupt\\CMIP5_breakpoints\\historical_rcp85qccceGW.txt")
      sys.argv.append("MIROC5r3i1p1")    
    #sys.argv.append( "MIROC5r1i1p1")# "ACCESS1-0r1i1p1")
    fn = sys.argv[1]
    #print all85.years()
    years=[]
    with open(fn, 'r') as f:
      lines = f.readlines()
      for line in lines:
        #print line
        (_,_,_,model,year,tas)=line.split()
        if not model in tasD:
          tasD[model]=[]
        tasD[model].append(float(tas))
        if years ==[]:
          modyr=model
          years.append(int(year))
        else:
          if model == modyr:
            years.append(int(year))
            
    ys=np.array([[[t,t],[t,t]] for t in tasD[sys.argv[2]]], dtype=np.float32)
    ys=tasD[sys.argv[2]]
 
    print ys
    breakpoints = {}
    realisations = {}
    for ii in range(1):
      xs=np.array([[[t,t],[t,t]] for t in [np.random.randn() for i in range(len(ys))]], np.float32)
      xs=np.array([i for i in range(len(ys))], np.float32)
      print np.shape(xs), np.shape(ys)
      bv=bivariate(ys[:], xs[:], critical = 12.5, anomalise=False) 
      a,b,c =bv.allPoints(debug=True, pr=0.01)
      print "run ", ii, "\nTi",a,"\nIx",b,"\nRange",c
      realisation = {}
      for row in range(len(a)-1):
        print "TiX", b[row]   
        key = str(b[row])
        try:
          if not key in breakpoints:
            breakpoints[key]=[]
            print "entered bp ", key
          breakpoints[key].append(str(a[row]))
        except Exception as e:
          print "exception ",e
          pass
    print breakpoints
    for y in np.sort(breakpoints.keys()):
      print y, len(breakpoints[y])
    sys.exit() 

#itewrative
    tasD={}
    fn = sys.argv[1]
    years=[]
    with open(fn, 'r') as f:
      lines = f.readlines()
      for line in lines:
        #print line
        (_,_,_,model,year,tas)=line.split()
        if not model in tasD:
          tasD[model]=[]
        tasD[model].append(float(tas))
        if years ==[]:
          modyr=model
          years.append(int(year))
        else:
          if model == modyr:
            years.append(int(year))
    ys=np.array(tasD[sys.argv[2]], dtype=np.float32)
    breakpoints = {}
    realisations = {}
    for i in range(2):
      xs=np.array([np.random.randn() for i in range(len(ys))], np.float32)
      bv=bivariate(ys[:], xs[:], critical = 12.5, anomalise=False) 
      a,b,c =bv.allPoints(debug=False)
      print "\nTi",a,"\nIx",b,"\nRange",c
      realisation = {}
      for row in range(len(a)-1):
        print b[row]   
        try:
          if not int(b[row]) in breakpoints:
            breakpoints[int(b[row])]=[]
          breakpoints[int(b[row])].append(float(a[row]))
        except:
          pass
    print breakpoints
    for y in np.sort(breakpoints.keys()):
      print y, years[y], len(breakpoints[y])
    sys.exit() 
#graphy bits
    tasD={}
    fn = sys.argv[1]
    years=[]
    with open(fn, 'r') as f:
      lines = f.readlines()
      for line in lines:
        #print line
        (_,_,_,model,year,tas)=line.split()
        if not model in tasD:
          tasD[model]=[]
        tasD[model].append(float(tas))
        if years ==[]:
          modyr=model
          years.append(int(year))
        else:
          if model == modyr:
            years.append(int(year))
    ys=np.array(tasD[sys.argv[2]], dtype=np.float32)
    xs=np.array([np.random.randn() for i in range(len(ys))], np.float32)  
    bv=bivariate(ys[:], xs[:], critical = 12.5, anomalise=False) 
    a,b,c =bv.allPoints(debug=False)
    print "\nTi",a,"\nIx",b,"\nRange",c
    
    ylist = {1900:0}
    
    for row in range(len(a)):
      try:
        print years[int(b[row])],a[row], b[row], c[0][row], c[1][row]
        ylist[years[int(b[row])+1]]= int(b[row])
      except:
        pass
    print np.shape(bv.Ti()), np.shape(years[:-1])
    ylist[2100]=len(ys)-1
    pline,=plt.plot(years[:-1],bv.Ti())
    
    plt.show()
    
    pl=[]
    ykeys=np.sort(ylist.keys())
    #print ykeys
    #print ylist
    fig, ax = plt.subplots()
    ax.set_title('%s %s' % (sys.argv[1], sys.argv[2]))
    pl.append(ax.plot(years, ys))
    for i in range(len(ykeys)-1):
      #print ykeys[i]
      #print ylist[ykeys[i]]
      sYs = np.array(ys[ylist[ykeys[i]]:ylist[ykeys[i+1]]+1])
      sXs = np.array(years[ylist[ykeys[i]]:ylist[ykeys[i+1]]+1])
      #print sXs, sYs
      stats=regress.analysed_regress(sYs,sXs)
      yhat, resid=regress.residuals(sYs, sXs, stats)
      pl.append(ax.plot(sXs, yhat, '-'))
      
    plt.show()
    sys.exit() 
    
    column=1
    for column in [1, 2, 3]:
        ys=np.array([d[column] for d in MyData], np.float32)    
        xs=np.array([np.random.randn() for i in range(len(MyData))], np.float32)   
        xs=np.array( [random.random()for d in MyData], np.float32)
        #xs = MyxData
        #print somedata.Header[i], somedata.Temps[i]
        bv=bivariate(ys[:], xs[:], critical = 12.5, anomalise=False) 
        #print somedata.Header[column],
        a,b,c =bv.allPoints(debug=False)
        print "\nTi",a,"\nIx",b,"\nRange",c
        for row in range(len(a)):
          try:
            print somedata.Header[column],a[row], b[row], c[0][row], c[1][row]
          except:
            pass
        #print bv.Ti()
        #print bv.nTi()
        
        #bv=bivariate(ys[0:36], xs[0:36], anomalise=False)
        #dumpme  = True
#       bv=bivariate(ys[86:9o], xs[86:90], anomalise=False) 
        if dumpme:
          print bv.maxIndexTi(), bv.maxTi()
          print "MEANNSES", 86, 90, bv.datamean(), bv.xmean(), bv.datastd(), bv.xstd()
          for k in range(4):
            try:
              print k+86, bv.ns()[k], bv.Sxy()[k], bv.Xi()[k],bv.Fi()[k], bv.Di()[k], bv.Ti()[k], bv.xi()[k], bv.yi()[k], bv.N()
            except:
              pass
    sys.exit() #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    column=1
    for column in [1,2,3]:
        ys=np.array([d[column] for d in MyData], np.float32)    
        xs=np.array([np.random.randn() for i in range(len(MyData))], np.float32)   
        xs=np.array( [random.random()for d in MyData], np.float32)
        xs = MyxData     
        #print somedata.Header[i], somedata.Temps[i]
        i = 0
        Ti = 100
        while Ti > 12.5:
            
            bv=bivariate(ys[i:], xs[i:], anomalise=False)        
            Ti=bv.maxTi()
            i=i+bv.maxIndexTi()
            print i, Ti, somedata.Header[column], MyData[i]
        
    column=1
    for column in [1,2,3]:
        ys=np.array([d[column] for d in MyData], np.float32)    
        xs=np.array([np.random.randn() for i in range(len(MyData))], np.float32)    
        #xs=np.array([random.random()for d in MyData], np.float32)    
        #print somedata.Header[i], somedata.Temps[i]
        i = len(MyData)             
        Ti = 100
        while Ti > 12.5:
            
            bv=bivariate(ys[:i], xs[:i], anomalise=False)        
            Ti=bv.maxTi()
            i=bv.maxIndexTi()
            print i, Ti, somedata.Header[column], MyData[i]
        
    print "Repeated with reinit"
    
    
    for column in [1,2,3]:
        ys=np.array([d[column] for d in MyData], np.float32)    
        xs=np.array([np.random.randn() for i in range(len(MyData))], np.float32)   
        bv=bivariate(ys[:], xs[:], anomalise=False)        
        Ti=bv.maxTi()
        i=bv.maxIndexTi()
        #print bv.Ti(), bv.maxTi()
        while Ti > 12.5:
            #print i
           # i = 4
           #bvv=bv.reinit(i, len(MyData)*np.ones(np.shape(bv.maxIndexTi()))) 
            bvv=bv.reinit(np.zeros(np.shape(bv.maxIndexTi())), i, 0.01, debug=False) 
            Ti = bvv.maxTi()
            i = bvv.maxIndexTi()
            #print bvv.nXi(), bvv.nYi(),bvv.mask(),bvv.nSxy(), bvv.nFi(),bvv.xi()[i], bvv.yi()[i], bvv.nns()[i], bvv.nNs(), bvv.nSxy()[i], bvv.nFi()[i], bvv.Fi()[-1]
            #sys.exit()
            #print i, Ti
            print i, Ti, somedata.Header[column], MyData[i]
    sys.exit()
    "/apollo/qccce/climate/qc1/cmip5/Downloads/ESG/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/mon/atmos/Amon/r1i1p1/tas/1/tas_Amon_CSIRO-Mk3-6-0_historical_r1i1p1_185001-200512.nc"
#    gwfn="../cmip5/tas_Amon_CSIRO-Mk3-6-0_historical_r1i1p1_185001-200512.GW"
#    fn="tas_Amon_CSIRO-Mk3-6-0_historical_r1i1p1_185001-200512.nc"
#    fn2="tas_Amon_CSIRO-Mk3-6-0_historical_r1i1p1_185001-200512_Ti.nc"
#    fn="tas_Amon_CSIRO-Mk3-6-0_piControl_r1i1p1_000101-050012.nc"
#    fn2="tas_Amon_CSIRO-Mk3-6-0_piControl_r1i1p1_000101-050012_Ti.nc"
#    fn="tas_Amon_CSIRO-Mk3-6-0_rcp85_r1i1p1_200601-210012.nc"
#    fn2="tas_Amon_CSIRO-Mk3-6-0_rcp85_r1i1p1_200601-210012_Ti.nc"
    gwfn="../cmip5/tas_Amon_CSIRO-Mk3-6-0_rcp85_r1i1p1_200601-210012.GW"
    fn="tas_Amon_CSIRO-Mk3-6-0_rcp85_r1i1p1_200601-210012.nc"
    fn2="tas_Amon_CSIRO-Mk3-6-0_rcp85_r1i1p1_200601-210012_Ti_GW.nc"


    nc=Nio.open_file(fn ,"r")
    shutil.copy(fn,fn2)
    tas=nc.variables["tas"][:]
    nc.close()
    xData = np.zeros((len(tas),))
    for i in range(len(xData)):
        xData[i] = np.random.randn()
    with open(gwfn, 'r') as gw:
        xd = []
        for line in gw.readlines():
            for i in range(12):
                xd.append(line.split()[1])
        xData=np.array(xd, np.float32)
    print xData 
    X=bivariate(tas, xData)
    nc=Nio.open_file(fn2 ,"rw")
    ncv=nc.create_variable("Ti", "f", nc.variables["tas"].dimensions)
    print np.shape(nc.variables["tas"]), np.shape(X.Ti())
    ncv[:]=np.array(X.Ti()[:], np.float32)
    ncv2=nc.create_variable("maxIndexTi", "i", nc.variables["tas"].dimensions[1:])
    
    ncv2[:]=np.array(X.maxIndexTi()[:], np.int32)[:]

    ncv3=nc.create_variable("maxTi", "f", nc.variables["tas"].dimensions[1:])
    ncv3[:]=np.array(X.maxTi()[:], np.float32)[:]
 
    nc.close()
