# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 14:40:29 2014

@author: s4493222
"""
import os
import numpy as np
import regress
from scipy.optimize import curve_fit
SVNRevision="$Revision: 307 $"
xdata = np.array([-2,-1.64,-1.33,-0.7,0,0.45,1.2,1.64,2.32,2.9])
ydata = np.array([0.699369,0.700462,0.695354,1.03905,1.97389,2.41143,1.91091,0.919576,-0.730975,-1.42001])

def func(x, p1,p2):
  return p1*np.cos(p2*x) + p2*np.sin(p1*x)

popt, pcov = curve_fit(func, xdata, ydata,p0=(1.0,0.2))

fn=os.environ["HOMEPATH"]+"\\Documents\\ReferenceData\\Rogers_Analysis_17Jul2014\\Data 4 Jim\\Global temps.csv"


def sigmoid(x, p1, p2, p3, p4):
  return p1 + p2/(1 + np.exp(p3 * (p4 - x)))

def sigmoid5(x, p1, p2, p3, p4, p5):
  return p1 + p2/(1 + np.exp(p3 * (p4 - x)))**p5

def sigmoid_est(x,y, five=False):
  try:
    slope, alpha=regress.regress(y, x) #get slope
    if slope > 0:
      p1 = np.min(y, axis=0)
      p2 = np.max(y, axis=0) - p1
      sign=1.0
    else:
      sign = -1.0
      p1 = np.max(y, axis=0)
      p2 = np.min(y, axis=0) - p1
    p4 = np.mean(x, axis=0)
    #print xprime, yprime,p2/(p1-yprime) -1,p4 - xprime, slope
    p3 =  sign * np.log(np.abs(slope))
  except Exception as e:
    print "Sig Est ", p1, p2, p3, p4
    raise e
  if five:
    #print "Est", (p1, p2, p3, p4, 1.0)
    return (p1, p2, p3, p4, 1.0)
  else:
    return (p1, p2, p3, p4)

def sig4(x,y):
  p = curve_fit(sigmoid, x, y, p0=sigmoid_est(x,y))[0]#,full_output=1)
  return sigmoid(x, p[0], p[1], p[2], p[3])

def sig5(x,y):
  p = curve_fit(sigmoid5, x, y, p0=sigmoid_est(x,y))[0]#,full_output=1)
  return sigmoid(x, p[0], p[1], p[2], p[3], p[4])

if __name__ == "__main__":
#call
  data=np.genfromtxt(fn,delimiter=",",names=True,filling_values =np.NaN) #fill with NaNs so must use their bounds
  print sigmoid_est(data["Year"][:-1], data["BCCA1r1"][:-1])
  p = curve_fit(sigmoid, data["Year"][:-1], data["BCCA1r1"][:-1], p0=sigmoid_est(data["Year"][:-1], data["BCCA1r1"][:-1],five=False))#,full_output=1)
  print p[0]
  params = p[0]
  
  
  import matplotlib.pyplot as plt
  fig = plt.figure()
  plt.plot(data["Year"][:-1], data["BCCA1r1"][:-1], 'b-')
  
  sig=sigmoid(data["Year"][:-1], params[0],params[1],params[2],params[3])
  plt.plot(data["Year"][:-1], sig, 'r-')

  plt.show()
  print sig4(data["Year"][:-1], data["BCCA1r1"][:-1])
