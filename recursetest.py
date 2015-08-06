# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 17:36:03 2014

@author: s4493222
"""

#!/usr/bin/env python
#support for the 
import numpy as np

import os
import sys
import bivariate_multi as bivariate
import lowess
from window40test_lowess import dataDict
import regress
import matplotlib.pyplot as plt
SVNRevision="$Revision: 297 $"

bcp85keys=['CCSM4r4i1p1', 'FGOALS-s2r2i1p1', 'HadGEM2-ESr4i1p1', 'IPSL-CM5A-LRr3i1p1', 'MIROC-ESM-CHEMr1i1p1', 'ACCESS1-0r1i1p1', 'FIO-ESMr2i1p1', 'GISS-E2-Rr1i1p1', 'FIO-ESMr1i1p1', 'GISS-E2-Rr1i1p3', 'bcc-csm1-1r1i1p1', 'MIROC5r2i1p1', 'EC-EARTHr9i1p1', 'inmcm4r1i1p1', 'CESM1-CAM5r2i1p1', 'FGOALS-s2r3i1p1', 'MIROC-ESMr1i1p1', 'MPI-ESM-LRr1i1p1', 'GISS-E2-Rr1i1p2', 'CSIRO-Mk3-6-0r3i1p1', 'CSIRO-Mk3-6-0r8i1p1', 'CNRM-CM5r1i1p1', 'CMCC-CMSr1i1p1', 'MRI-CGCM3r1i1p1', 'CSIRO-Mk3-6-0r6i1p1', 'MPI-ESM-LRr2i1p1', 'GISS-E2-Hr1i1p2', 'CESM1-BGCr1i1p1', 'HadGEM2-AOr1i1p1', 'GISS-E2-Hr1i1p1', 'EC-EARTHr12i1p1', 'GFDL-CM3r1i1p1', 'GISS-E2-Hr1i1p3', 'EC-EARTHr8i1p1', 'GFDL-ESM2Mr1i1p1', 'CSIRO-Mk3-6-0r1i1p1', 'HadGEM2-CCr1i1p1', 'EC-EARTHr1i1p1', 'MPI-ESM-LRr3i1p1', 'CCSM4r5i1p1', 'CanESM2r5i1p1', 'HadGEM2-ESr1i1p1', 'HadGEM2-ESr2i1p1', 'FGOALS-s2r1i1p1', 'GFDL-ESM2Gr1i1p1', 'CCSM4r6i1p1', 'ACCESS1-3r1i1p1', 'CNRM-CM5r2i1p1', 'MIROC5r1i1p1', 'IPSL-CM5B-LRr1i1p1', 'CSIRO-Mk3-6-0r2i1p1', 'FIO-ESMr3i1p1', 'MIROC5r3i1p1', 'CCSM4r3i1p1', 'CMCC-CESMr1i1p1', 'CNRM-CM5r4i1p1', 'CSIRO-Mk3-6-0r4i1p1', 'CSIRO-Mk3-6-0r9i1p1', 'NorESM1-Mr1i1p1', 'EC-EARTHr2i1p1', 'CESM1-CAM5r3i1p1', 'FGOALS-g2r1i1p1', 'CCSM4r1i1p1', 'bcc-csm1-1-mr1i1p1', 'CNRM-CM5r6i1p1', 'IPSL-CM5A-LRr1i1p1', 'CanESM2r1i1p1', 'CSIRO-Mk3-6-0r5i1p1', 'IPSL-CM5A-MRr1i1p1', 'MPI-ESM-MRr1i1p1', 'CNRM-CM5r10i1p1', 'CESM1-CAM5r1i1p1', 'CSIRO-Mk3-6-0r7i1p1', 'CSIRO-Mk3-6-0r10i1p1', 'CanESM2r4i1p1', 'CanESM2r3i1p1', 'HadGEM2-ESr3i1p1', 'IPSL-CM5A-LRr4i1p1', 'CanESM2r2i1p1', 'IPSL-CM5A-LRr2i1p1', 'CCSM4r2i1p1', 'NorESM1-MEr1i1p1', 'CMCC-CMr1i1p1', 'BNU-ESMr1i1p1']

class recurseTestException(Exception):
  def __str__(self):
    return "recurseTestException "+self.message 
    
class recurse(object):
  def __init__(self, ys, xs, years, model, pr=0.01, smooth=False, anom=False, onethreshold=True, trim=1, debug=False, ConstSxy=True, withshifts=False):
    try:
      #print "Pr size", pr, size
      self.__size=len(xs)
      #print self.__size, type(self.__size)
      self.__threshold=bivariate.critTi(pr, self.__size)
      self.__breakpoints=np.zeros((self.__size,))
      self.__breakyears = {}
    except Exception as e:
      print bivariate.ExceptionInfo()
      print "__init__",str(e)
      raise recurseTestException("__init__:"+str(e))
      
    if smooth:
      if anom:
        #print "window smooth anom"
        txs, tys, td1, td2 = lowess.R().lowess(np.array(xs),np.array(ys), f=1./4., iter=1)
        self.__ys = np.array(ys)-np.array(tys)      
        self.__xs=np.copy(xs)
      else:
        #print "window smooth not anom"
        txs, tys, td1, td2 = lowess.R().lowess(np.array(xs),np.array(ys), f=1./4., iter=1)
        self.__ys = np.array(ys)              
        self.__xs=np.array(tys)+np.array(xs)
    else:
      #print "window not smooth"
      self.__ys=np.array(ys)
      self.__xs=np.copy(xs)
    #print "window call bivariate"
    try:
      #self.__bv=bivariate.bivariate(self.__ys, self.__xs, critical=self.__threshold, anomalise=False,constantsxy=ConstSxy)
      self.__bv=bivariate.bivariate(self.__ys, self.__xs, critical=None, pr=pr, anomalise=False,constantsxy=ConstSxy)
      if withshifts:
        ap_MaxTis, ap_MaxIndexes, (ap_lows, ap_highs), ap_shifts = self.__bv.allPoints(pr, withshifts=True)
      else:
        ap_MaxTis, ap_MaxIndexes, (ap_lows, ap_highs)= self.__bv.allPoints(pr, withshifts=False)
      #print "BV.BIVARIATE -> ", len(ap_MaxTis), len(ap_MaxIndexes), len(ap_lows), len(ap_highs)
    except Exception as e:
      print bivariate.ExceptionInfo()
      raise recurseTestException("recurse.__init__ call to bivariate:"+str(e))  
    #print ap_MaxTis, ap_MaxIndexes, (ap_lows, ap_highs)
    
    
    if trim == 1:  
#==============================================================================
##       This mode implements a trimming mode whereby breakpoints tested by
##    bracketing them and accepting them if there ids a breakpoint in the interval
##    The alternative (mode 2) is to replace them with the new breakpoint    
#==============================================================================
    
      try:
        stepn = 0
        if debug: print "trimming"
        trimmed = True
        while trimmed:
          trimmed = False
          for mi in range(len(ap_MaxIndexes)):
            if debug: print "considering mi"
            if not trimmed:
              stepn = 1
              self.__bv.reinit(ap_lows[int(mi)], ap_highs[int(mi)]+1, pr)
              stepn=2
              #print self.__bv.maxnTi() , ap_lows[int(mi)], ap_MaxIndexes[int(mi)], ap_highs[int(mi)], bivariate.critTi(pr, 1 + ap_highs[int(mi)]-ap_lows[int(mi)])
              stepn = 3
              if self.__bv.maxnTi() < bivariate.critTi(pr, 1 + ap_highs[int(mi)]-ap_lows[int(mi)]):
                stepn = 4
                trimmed = True
                if debug: print "trimmer removed ", mi, ap_MaxIndexes[int(mi)], " between ", ap_lows[int(mi)], ap_highs[int(mi)]
                stepn = 5
                #print mi, len(ap_MaxTis), len(ap_MaxIndexes), len(ap_lows), len(ap_highs)
                #print ap_MaxTis, " becomes ", 
                pop1 = ap_MaxTis.pop(int(mi))
                #print ap_MaxTis
                pop2 = ap_MaxIndexes.pop(int(mi))
                pop3 = ap_lows.pop(int(mi))
                pop4 = ap_highs.pop(int(mi))
                with open("removals.txt","a") as rmf:
                  rmf.write("from %s trimmer removed element %s (Ti:%s ) %s between %s and %s crit=%s \n" % (str(model), str(mi), str(pop1), str(pop2), str(pop3), str(pop4), str( bivariate.critTi(pr, 1 + ap_highs[min(len(ap_highs)-1,mi)]-ap_lows[min(len(ap_lows)-1,mi)]))))
      except Exception as e:
        raise recurseTestException("at trim state "+str(stepn)+" in __init__:"+str(e))

    if trim == 2:  
#==============================================================================
##       This mode implements a trimming mode whereby breakpoints tested by
##    bracketing them and accepting them if there ids a breakpoint in the interval
##    The alternative (mode 2) is to replace them with the new breakpoint    
#==============================================================================
    
      try:
        stepn = 0
        if debug: print "trimming mode 2"
        trimmed = True
        while trimmed:
          trimmed = False
          for mi in range(len(ap_MaxIndexes)):
            if debug: print "considering mi"
            if not trimmed:
              stepn = 1
              self.__bv.reinit(ap_lows[int(mi)], ap_highs[int(mi)]+1, pr)
              stepn=2
              print self.__bv.maxnTi() , ap_lows[int(mi)], ap_MaxIndexes[int(mi)], ap_highs[int(mi)], bivariate.critTi(pr, 1 + ap_highs[int(mi)]-ap_lows[int(mi)])
              if self.__bv.maxnTi() < bivariate.critTi(pr, 1 + ap_highs[int(mi)]-ap_lows[int(mi)]):  
                trimmed = True
                if debug: print "trimmer (mode 2) removed ", mi, ap_MaxIndexes[int(mi)], " between ", ap_lows[int(mi)], ap_highs[int(mi)]
                pop1 = ap_MaxTis.pop(int(mi))
                pop2 = ap_MaxIndexes.pop(int(mi))
                pop3 = ap_lows.pop(int(mi))
                pop4 = ap_highs.pop(int(mi))
                with open("removals.txt","a") as rmf:
                  rmf.write("from %s trimmer (mode 2) removed element %s (Ti:%s ) %s between %s and %s crit=%s \n" % (str(model), str(mi), str(pop1), str(pop2), str(pop3), str(pop4), str( bivariate.critTi(pr, 1 + ap_highs[min(len(ap_highs)-1,mi)]-ap_lows[min(len(ap_lows)-1,mi)]))))
              elif self.__bv.maxIndexnTi() != ap_MaxIndexes[int(mi)]:
                trimmed = True
                pop1 = ap_MaxTis[int(mi)]
                pop2 = ap_MaxIndexes[int(mi)]
                pop3 = ap_lows[int(mi)]
                pop4 = ap_highs[int(mi)]
                ap_MaxTis[int(mi)] = self.__bv.maxnTi()               
                ap_MaxIndexes[int(mi)] = self.__bv.maxIndexnTi()
                with open("removals.txt","a") as rmf:
                  rmf.write("from %s trimmer (mode 2) substituted element %s (Ti:%s ) %s with %s between %s and %s crit=%s \n" % (str(model), str(mi), str(ap_MaxIndexes[int(mi)]), str(pop1), str(pop2), str(pop3), str(pop4), str( bivariate.critTi(pr, 1 + ap_highs[int(mi)]-ap_lows[int(mi)]))))
                
      except Exception as e:
        raise recurseTestException("at mode 2 trim state "+str(stepn)+" in __init__:"+str(e))

    try:  
      stepn = 0
      for mi in range(len(ap_MaxIndexes)):
        #print mi
        stepn = 1
        if ((onethreshold and (ap_MaxTis[int(mi)] >= self.__threshold)) or 
            (not onethreshold and (ap_MaxTis[int(mi)] >= 
                                  bivariate.critTi(pr, 1 + ap_highs[int(mi)]-ap_lows[int(mi)])))):
          self.__breakpoints[int(ap_MaxIndexes[int(mi)])] += 1
          stepn = 2
          self.__breakyears[years[int(ap_MaxIndexes[int(mi)])]] = None
          stepn = 3
          if withshifts:
            self.__breakyears[years[int(ap_MaxIndexes[mi])]] = (ap_MaxTis[mi], years[int(ap_lows[mi])], years[int(ap_highs[mi])-1], ap_shifts[mi])
          else:
            self.__breakyears[years[int(ap_MaxIndexes[mi])]] = (ap_MaxTis[mi], years[int(ap_lows[mi])], years[int(ap_highs[mi])-1])
          
    except Exception as e:
      print "Exception -----------------------------------------------------------"
      try:
        print "mi",mi
        #print "years",years
        print "ap_MaxIndexes",ap_MaxIndexes
        print "ap_MaxTis",ap_MaxTis
        print "ap_lows",ap_lows
        print "ap_highs",ap_highs
        print "self.__breakyears",self.__breakyears
        print "years[int(ap_MaxIndexes[int(mi)])]",years[int(ap_MaxIndexes[int(mi)])]
      except:
        pass
      print bivariate.ExceptionInfo()
      raise recurseTestException("Reporting loop of __init__: "+str(stepn)+" "+str(e))
      
      
  def breakpoints(self):
    return self.__breakpoints

  def breakyears(self):
    return self.__breakyears

  def plot(self, x, y, years, title=0.0, save=True):
    ys = self.__ys
    xs = self.__xs
    print x, y, years
    #print self.__breakyears.keys()
    breaks=[years[0]]
    breaks.extend(np.sort(self.__breakyears.keys()))
    breaks.extend([years[-1]])
    #print breaks
    fig, ax = plt.subplots()
    if title != 0.0:
      ax.set_title(sys.argv[2]+ " "+title)
    pl = []
    #print "LOWESS",lowess(xs,ys)
    #sys.exit()
    #print "HERE"
    pl.append(ax.plot(years, y))

    for i in range(len(breaks)-1):
      try:
        pt = 0
        sYs = np.array(y[breaks[i]-breaks[0]:breaks[i+1]+1-breaks[0]])
        sXs = np.array(years[breaks[i]-breaks[0]:breaks[i+1]+1-breaks[0]])
        #print breaks[i]-breaks[0], breaks[i+1]+1-breaks[0], sXs, sYs
        #print "JIM", len(xs), len(ys)
        pt = 1
        stats=regress.analysed_regress(sYs,sXs)
        pt = 2
        yhat, resid=regress.residuals(sYs, sXs, stats)
        pt = 3
        pl.append(ax.plot(sXs, yhat, '-'))
      except Exception as e:
        print "Exception at "+str(e),pt
        
    if not save:
      plt.show()
    else:
      plt.savefig(sys.argv[2]+ " "+title+'.png')
    return breaks

    #JHR 7-5-2014 support for selection of breakpoints given an accelerating series, to avoid over reaching
  #def validslice(self, MaxIndexes):
       
def process(dataD, gcm, filename, smooth=True, anom=True):
    ys=np.array(dataD[gcm], dtype=np.float32)
    if os.path.exists(filename):
      os.remove(filename)
    for j in range(100):
      xs=np.array([np.random.randn() for i in range(len(ys))], np.float32)
      #print "process window"
      bp=recurse(ys, xs, gcm, smooth=smooth, anom=anom).breakpoints()
      #print "process segment"
      with open(filename, "a") as outf:
        for line in bp:
          outf.write(str(line)+"\n")
     
			      
if __name__ == "__main__":
#itewrative
    tasD={}
    if len(sys.argv) <=1:
      sys.argv.append(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\CMIP5_breakpoints\\historical_rcp85qccceGW.txt")
    #sys.argv.append( "MIROC5r1i1p1")# "ACCESS1-0r1i1p1")
    fn = sys.argv[1]
    all85=dataDict(fn)
    #print all85.years()
    sys.argv.append("MIROC5r3i1p1")    
    if len(sys.argv) <=2:
      for model in all85.models():
        if not os.path.exists(model+'.summary.txt'):
          sys.argv.append(model)
      #sys.argv.append( "MIROC5r3i1p1")# "ACCESS1-0r1i1p1")
      #sys.argv.append( "MIROC5r2i1p1")# "ACCESS1-0r1i1p1")
      #sys.argv.append( "MIROC5r1i1p1")# "ACCESS1-0r1i1p1")
      #sys.argv.append( "ACCESS1-0r1i1p1")
    years = all85.years()[sys.argv[2]]
    print "YEARS",sys.argv[1], sys.argv[2], years
    tasD = all85[sys.argv[2]]
    model = sys.argv[2]
    ys=np.array(tasD, dtype=np.float32)
    byssummary = []
    if os.path.exists(model+'.years.txt'): os.remove(model+'.years.txt')
    if os.path.exists(model+'.bp.txt'): os.remove(model+'.bp.txt')
    for j in range(100):  ###############################################################
      xs=np.array([np.random.randn() for i in range(len(ys))], np.float32)
      
      try:
        stepn=0
        bpv=recurse(ys, xs, years, "%s_run_%d" % (model, j), 0.01, smooth=False, anom=False, debug=True, trim=1)
        bp = bpv.breakpoints()
        bys = bpv.breakyears()
        stepn=1
        byfn=model+'.years.txt'
        bpfn=model+'.bp.txt'     
        with open(byfn, "a") as byf:
          for k in np.sort(bys.keys()):
            print k, bys[k]
            byf.write("%s %s\n" % (str(k), str(bys[k])))
          
        with open(bpfn, "a") as bpf:
          for p in bp:
            bpf.write("%s\n" % (str(p),))
        
        byssummary.append(bpv.plot(xs, ys, years, str(j)))
      except Exception as e:
        print "MAIN Exception Step ", stepn, str(e)
        pass
      
    #print byssummary  
    byssummarycounts={} #we will count how many actual different predictions we get
    
    plt.show()
    for i in range(len(byssummary)):
      ax1=plt.plot(byssummary[i], np.ones(np.shape(byssummary[i]))*i, "bo")
      bskey=str(byssummary[i])
      if not (bskey in byssummarycounts):
        byssummarycounts[bskey] = 0
      byssummarycounts[bskey] += 1
    #plot the variance, skew and kurtosis all normalised between 0 and 100
    import runningstats
    s=runningstats.runningMoments(ys)
    sdelta=runningstats.runningMoments(ys[1:]-ys[:-1])
      
    (mean, variance, sigma, skew, kurtosis) = s.moments()
    #print (mean, variance, sigma, skew, kurtosis)
    vararray= [variance]
    skewarray= [skew]
    kurtarray= [ kurtosis]
    while s.step():
      (mean, variance, sigma, skew, kurtosis) = s.moments()
      #print (mean, variance, sigma, skew, kurtosis)      
      vararray.append(variance)
      skewarray.append(skew)
      kurtarray.append(kurtosis)

    vararray = np.array(vararray)
    skewarray = np.array(skewarray)
    kurtarray = np.array(kurtarray)
    vararray = 50 + 50  * (np.array(vararray)-np.mean(vararray))/(np.max(vararray)-np.min(vararray))
    skewarray = 50 + 50.* (np.array(skewarray)-np.mean(skewarray))/(np.max(skewarray)-np.min(skewarray))
    kurtarray = 50 + 50.* (np.array(kurtarray)-np.mean(kurtarray))/(np.max(kurtarray)-np.min(kurtarray))
    ax3=plt.plot(years[15:15+len(vararray)], vararray,"g--")
    ax4=plt.plot(years[15:15+len(vararray)], skewarray,"r--")
    ax5=plt.plot(years[15:15+len(vararray)], kurtarray,"k--")
    
    (mean, variance, sigma, skew, kurtosis) = sdelta.moments()
    #print (mean, variance, sigma, skew, kurtosis)
    vararray= [variance]
    skewarray= [skew]
    kurtarray= [ kurtosis]
    while sdelta.step():
      (mean, variance, sigma, skew, kurtosis) = sdelta.moments()
      #print (mean, variance, sigma, skew, kurtosis)      
      vararray.append(variance)
      skewarray.append(skew)
      kurtarray.append(kurtosis)

    vararray = np.array(vararray)
    skewarray = np.array(skewarray)
    kurtarray = np.array(kurtarray)
    vararray = 50 + 50  * (np.array(vararray)-np.mean(vararray))/(np.max(vararray)-np.min(vararray))
    skewarray = 50 + 50.* (np.array(skewarray)-np.mean(skewarray))/(np.max(skewarray)-np.min(skewarray))
    kurtarray = 50 + 50.* (np.array(kurtarray)-np.mean(kurtarray))/(np.max(kurtarray)-np.min(kurtarray))
    ax3=plt.plot(years[15:15+len(vararray)], vararray,"g-")
    ax4=plt.plot(years[15:15+len(vararray)], skewarray,"r-")
    ax5=plt.plot(years[15:15+len(vararray)], kurtarray,"k-")
    
    rescaled = np.array(100.0*(ys-np.min(ys))/(np.max(ys)-np.min(ys)))
    ax6=plt.plot(years, rescaled,"b-")
    plt.title= sys.argv[2]
    plt.savefig(sys.argv[2]+ "_Summary.png")
    maxcount = 0
    maxlist = []
    
    with open(model+'.summary.txt', "w") as bysf:      
      for k in np.sort(byssummarycounts.keys()):
        if byssummarycounts[k] == maxcount:
          maxlist.append(k)
        elif byssummarycounts[k] > maxcount:
          maxlist = [k]
          maxcount = byssummarycounts[k]
        print k, byssummarycounts[k]
        bysf.write("%s %s\n" % (str(k), str(byssummarycounts[k])))
    pl = []
    for breakk in maxlist:
      breaks=[int(b) for b in breakk[1:-1].split(',')]
      #print "BREAKS",breaks
      for i in range(1,len(breaks)):
        try:
          pt = 0
          #print "RANGE",years[0],breaks[i-1],years[0],breaks[i], range(years[0]-breaks[i-1],years[0]-breaks[i])
          sYs = np.array(rescaled[breaks[i-1]-years[0]:breaks[i] - years[0]])
          sXs = np.array(range(breaks[i-1],breaks[i]))
         # print "SXY",sYs, sXs
          #print breaks[i]-breaks[0], breaks[i+1]+1-breaks[0], sXs, sYs
          #print "JIM", len(xs), len(ys)
          pt = 1
          stats=regress.analysed_regress(sYs,sXs)
          pt = 2
          yhat, resid=regress.residuals(sYs, sXs, stats)
          pt = 3
          pl.append(plt.plot(sXs, yhat, '-'))
        except Exception as e:
          print bivariate.ExceptionInfo()
          print "Exception at ",str(e),pt
      
