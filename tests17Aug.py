# -*- coding: utf-8 -*-
"""
Created on Sat Jul 19 12:32:03 2014

@author: s4493222
"""
import os
import numpy as np
import bivariate_multi as bivariate
import recursetest
import random
#import lowess
import datetime
import matplotlib.pylab as plt
import regress
#import math
import AIC
import copy
import scipy.stats as stats
SVNRevision="$Revision: 286 $"
fn=os.environ["HOMEPATH"]+"\\Documents\\ReferenceData\\Rogers_Analysis_17Jul2014\\GISSTEMPV3.csv"
#fn=os.environ["HOMEPATH"]+"\\Documents\\ReferenceData\\Rogers_Analysis_17Jul2014\\Data 4 Jim\\Global temps.csv"
fn=os.environ["HOMEPATH"]+"\\Documents\\abrupt\\CMIP5_breakpoints\\Artyfishul_dartar4Jim (1).csv"
#fn=os.environ["HOMEPATH"]+"\\Documents\\abrupt\\CMIP5_breakpoints\\Artyfishul_dartar4Jim.csv"
TraceFile = "genericTrace.txt" #this name is replace at run time.

class RCTestException(Exception):
  def __init__(self,msg):
    self.msg=("RC Test"+str(msg))
  def __str__(self):
    return str(self)  

trimmode=0

rogerslist={
'Year':[],
'NHem':[1925,1987, 1997],
'24S24N':[1937,1972,1997], 
'Glob':[1929, 1979, 1996], 
'SHem':[1936, 1969, 1996], 
'90S24S':[], 
'24SEQU':[], 
'44N64N':[1920,1988,1997], 
'24N90N':[], 
'44S24S':[1937, 1972, 1997], 
'EQU24N':[1926, 1987], 
'64S44S':[1937, 1968], 
'24N44N':[1937, 1964, 1987, 1998], 
'90S64S':[1955], 
'64N90N':[1920, 1955, 1988, 2005]
}

ConstSxy = True #Introduced to see whether Roger's version of the bivariate test was different in early runs. Decision:Use as is.
#jims 0 was withoffset=false, scoffset = 0 plus
#jims 2 was inconsistently coded
#jims 3 was withoffset=true, scoffset = 0
#jims 4 was withoffset=true, scoffset = 1  
RunName="JimsTest_7"
withoffset=True  # Initially I was shooting one year shorter in the reduction. Set to True to fix this.Tested and True is it.
#Now see if strucchnage is reporting year of chnage of first year chnaged
scoffset = 0 #set to zero for initially coded case - conclusion - 0 better than -1 and +1
runlen=100

if ConstSxy:
  RunName += "_ConstantSxy"
else:
  RunName += "_ProgressiveSxy"
 
    
def yearlist(Years, Counts):
  l = []
  for x in range(len(Counts)):
    if Counts[x] != 0.0:
      l.append(Years[x])
  return l
  
def valid(vect):
 # print vect
  lo=0
  hi=len(vect)-1
  while lo < hi and np.isnan(vect[lo]):
#    print "+",Years[lo], vect[lo], Years[hi], vect[hi]
    lo +=1
  while hi >=lo and np.isnan(vect[hi]):
#    print "i",Years[lo], vect[lo], Years[hi], vect[hi]
    hi -= 1
  #print Years[lo], vect[lo], Years[hi], vect[hi]
  return lo, hi+1
  

def plotbreaks(x, y, breakpoints, title=None, subtitle="", breakyears=None, offset=0):
  fig = plt.figure()
  breaklist = []
  statslist = []
  for i in range(len(breakpoints)):
    if breakpoints[i] != 0:
      breaklist.append(i+offset+1)
  plt.plot(x, y, 'b-')
  lo = 0
  hi = len(y)
  #now perform and plot regressions
  for b in breaklist:
    print b
    stats = regress.analysed_regress(y[lo:b], x[lo:b])
    yh, r = regress.residuals(y[lo:b], x[lo:b], stats)
    plt.plot(x[lo:b], yh, 'r-')
    statslist.append(stats)
    lo = b
  stats = regress.analysed_regress(y[lo:hi], x[lo:hi])
  statslist.append(stats)
  yh, r = regress.residuals(y[lo:hi], x[lo:hi], stats)
  plt.plot(x[lo:hi], yh, 'r-')
  
  #plot the list of years also (this is probably from the strucchange analysis)
  if breakyears != None:
    lo = 0
    hi = len(y)
    breaklist2 = []
    statslist2 = []
    for yr in breakyears:
      #x is years
      breaklist2.append(list(x[:]).index(yr)+scoffset)
    #now perform and plot regressions
    for b in breaklist2:
      print title, "list 2", b
      stats = regress.analysed_regress(y[lo:b], x[lo:b])
      statslist2.append(stats)
      yh, r = regress.residuals(y[lo:b], x[lo:b], stats)
      plt.plot(x[lo:b], yh, 'g--')
      lo = b
    stats = regress.analysed_regress(y[lo:hi], x[lo:hi])
    statslist2.append(stats)
    yh, r = regress.residuals(y[lo:hi], x[lo:hi], stats)
    plt.plot(x[lo:hi], yh, 'g--')
#plot roger's breaklist as well
    
    
  breaklist3=[]
  statslist3=[]
  for b in rogerslist[title]:
    breaklist3.append(list(x[:]).index(b))
  lo = 0
  hi = len(y)
  for b in breaklist3:
    print title,"list 3",b
    stats = regress.analysed_regress(y[lo:b], x[lo:b])
    statslist3.append(stats)
    yh, r = regress.residuals(y[lo:b], x[lo:b], stats)
    plt.plot(x[lo:b], yh, 'r--')
    lo = b
  stats = regress.analysed_regress(y[lo:hi], x[lo:hi])
  statslist3.append(stats)
  yh, r = regress.residuals(y[lo:hi], x[lo:hi], stats)
  plt.plot(x[lo:hi], yh, 'r--')
    
  

  fig.suptitle(RunName+" "+title+"\n"+subtitle)
  fig.savefig(RunName+" "+title+'.png')  
  plt.show()
  return statslist, statslist2, statslist3

def plotbreaksSummary(x, y, breaks, title=None, subtitle="", breakyears=None, breakyearsAve=None, offset=0):
  fig = plt.figure()
  breaklist = []
  statslist = []
  for yr in breaks:
    #x is years
    breaklist.append(list(x[:]).index(yr)+scoffset+1)
  #now perform and plot regressions

  plt.plot(x, y, 'b-')
  lo = 0
  hi = len(y)
  #now perform and plot regressions
  
  if breakyearsAve != None:
    lo = 0
    hi = len(y)
    breaklist4 = []
    statslist4 = []
    for yr in breakyearsAve:
      #x is years
      breaklist4.append(list(x[:]).index(yr)+scoffset)
    #now perform and plot regressions
    for b in breaklist4:
      print title, "list 2", b
      stats = regress.analysed_regress(y[lo:b], x[lo:b])
      statslist4.append(stats)
      yh, r = regress.residuals(y[lo:b], x[lo:b], stats)
      plt.plot(x[lo:b], yh, 'k--')
      lo = b
    stats = regress.analysed_regress(y[lo:hi], x[lo:hi])
    statslist4.append(stats)
    yh, r = regress.residuals(y[lo:hi], x[lo:hi], stats)
    plt.plot(x[lo:hi], yh, 'k--')
#plot roger's breaklist as well
    
  if breakyears != []:  
    breaklist3=[]
    statslist3=[]
    lo = 0
    hi = len(y)
    for yr in breakyears:
      #x is years
      breaklist3.append(list(x[:]).index(yr)+scoffset)
    if breaklist3 != []:
      for b in breaklist3:
        print title,"list 3",b
        stats = regress.analysed_regress(y[lo:b], x[lo:b])
        statslist3.append(stats)
        yh, r = regress.residuals(y[lo:b], x[lo:b], stats)
        plt.plot(x[lo:b], yh, 'r--')
        lo = b
      stats = regress.analysed_regress(y[lo:hi], x[lo:hi])
      statslist3.append(stats)
      yh, r = regress.residuals(y[lo:hi], x[lo:hi], stats)
      plt.plot(x[lo:hi], yh, 'r--')
    
  

  fig.suptitle(RunName+" "+title+"\n"+subtitle)
  fig.savefig(RunName+" "+title+'.png')  
  plt.show()
  return statslist, statslist3

def reposition_break(testdata, controldata ,datayears, breakyear, model):
  rb = datayears.index(breakyear)
  #now do three regressions. all of sample and left and right of break
  
def resample_break(testdata, datayears, N=30):
  brks=[]
  tis=[]
  for i in range(N):
    controldata = np.array([random.random() for y in datayears])
    bv=bivariate.bivariate(testdata,controldata, anomalise=False, pr=0.01)
    tis.append(bv.maxTi())
    brks.append(datayears[bv.maxIndexTi()])
  return stats.norm.fit(brks), stats.norm.fit(tis)   

def convergentBreaks_Inner(testdata, controldata, datayears, aicControl, model, trace=True, shallow=False, keepFirst=False):
  '''
  This code attempst to iteratively test all 
  '''
  #The issue with convergence is that whether a breakpoint is admitted to the yearly breaks is determined 
  #initial analysis over full data
  if trace:
      tf=open(TraceFile,"a")

  rt0=recursetest.recurse(testdata, controldata, datayears, model, smooth=False, trim = 0, pr=0.01, anom=False)
  
  oyears = rt0.breakyears() #initial set of breakpoints 
  byears = [datayears[0]]
  byears.extend([k for k in np.sort(oyears.keys())])
  #print "BYEARS", byears
  #byears=np.insert(byears, 0, datayears[0])
  
  if trace:
    print >>tf,model, "trace=",trace, "shallow=",shallow
    for i in range(len(testdata)):
      print >>tf,i, testdata[i], controldata[i], datayears[i], aicControl[i]
  byears.append(datayears[-1]) 
  if trace: print >>tf,"BYEARS", byears
  initialBreaks= copy.copy(byears)
  low=0
  crit30 = bivariate.critTi(0.01, 30)
  print "initially", np.sort(oyears.keys()), crit30
  testedlist = []
  statlist = []
  newbreaks=[datayears[0]] #the first list start is preserved
  fails = 0
  while len(byears) > 2:
    popped=byears.pop(0)
    testedlist.append(popped)
    lo=np.argwhere(datayears==newbreaks[-1])[0][0]
    hi=np.argwhere(datayears==byears[1])[0][0]+1
    testyr = byears[0]    
    print "try", datayears[lo], datayears[hi-1],
    rt1=recursetest.recurse(testdata[lo:hi], controldata[lo:hi], datayears[lo:hi], model, smooth=False, trim = 0, pr=0.01, anom=False)
    print rt1.breakyears().keys()
    nlist=rt1.breakyears().keys()
    nlist.sort()
    if fails > 0:
      print nlist, "shortened to ",nlist[fails:]
    nlist = nlist[fails:]
    lennlist=len(nlist)
    nlist.extend(byears[1:])
    
    byears=nlist[:]
    if lennlist == 0:
      print popped, " gone, revised list empty"
    else:
      if lennlist > 1: #then lo and hi need to be computed
        hi=np.argwhere(datayears==byears[1])[0][0]+1       
        testyr =byears[0] 
      ystats, tstats = resample_break(testdata[lo:hi], datayears[lo:hi])
      print popped, datayears[hi-1],ystats, tstats, byears[0],
      if tstats[0] + 2 * tstats[1] < crit30:
        print " **"
        if lennlist > 0:
          fails += 1
          
      else:
        if crit30 > tstats[0]:
          print " * "
        else:
          print "   ",
        if abs(testyr - ystats[0] ) -0.5 <= 2 * ystats[1]:
          print "YOK",
          #KeepOldDropNew
          newbreaks.append(round(ystats[0],0))
          statlist.append((ystats, tstats))
          fails = 0
        else:
          print " Y* ",
          #So does the spread of possible break years actually exceed the possible bounds - this is a sign of instability or trend, or badness, or we may already have this point eqarlier
          if ystats[0] - 2 * ystats[1] < min(testyr, datayears[lo]) or ystats[0] + 2 * ystats[1] > max(testyr, datayears[hi-1]):
            print" drop", testyr,"likehotpotato"
            if lennlist > 0: fails += 1
          else:
            byears[0] = round(ystats[0],0)
            newbreaks.append(byears[0])
            statlist.append((ystats, tstats))
            print "revised",testyr,"to",byears[0]
            fails = 0
  newbreaks.append(datayears[-1])
  
  print "\nReturning",initialBreaks, "->", newbreaks, statlist
  return initialBreaks, newbreaks, statlist
  
def collate(discoveryList, cutoff = 0.9, reps=100):
#this code collects together the discovered years         
  allResults = {}
  for model in discoveryList:
    result = {}
    for br in discoveryList[model]:
    #now process the sets of years and produce stats
      N=discoveryList[model][br][0]
      for yr in br:
        if not yr in result:
          result[yr] = float(N)
        else:
          result[yr] +=N
    for yr in result:
      result[yr] /= float(reps)
   #result now has total numbers found for each year   

    keys = list(np.sort(result.keys()))
    if keys != []:
      k1 = keys[0]
      count = 1
      tot = result[k1]
      ktot = k1 * tot
      
      final = []  
      for k in keys[1:]:
        if k-k1 <=1:
          val = result[k]
          tot += val
          ktot += k * val
          count +=1
        else:
          if count > 0:
            ave=tot
            year = ktot/tot
            if ave >= cutoff:
              final.append((year, ave))
          count = 1
          tot = result[k]
          ktot = k * tot
        k1 = k
      if count !=0.:
          ave=tot
          year = ktot/tot
          count = 0
          ktot = 0.
          tot = 0.
          if ave >= cutoff:
            final.append((year, ave))
    else:
      final= []
    #print final  
    allResults[model] = final
  return allResults

controlmodes=set(["control","years","flat","trend"])
def convergentBreaks(testdata, controldata, datayears, model, mode="years", guide="AIC"):
  if not mode in controlmodes:
    raise RCTestException("keyword 'mode' should be in "+str(controlmodes))
  if guide == "AIC":
    if mode=="flat":
      ctrl=np.ones(np.shape(controldata))
      return convergentBreaks_Inner(testdata, controldata, datayears, ctrl, model)
    elif mode=="years":
      return convergentBreaks_Inner(testdata, controldata, datayears, datayears, model)
    elif mode=="control":
      return convergentBreaks_Inner(testdata, controldata, datayears, controldata, model)
    else:
      stats = regress.analysed_regress(testdata, controldata)
      yh, r = regress.residuals(testdata, controldata, stats)
      return convergentBreaks_Inner(testdata, controldata, datayears, yh, model)
  
if __name__ == "__main__":
 
  data=np.genfromtxt(fn,delimiter=",",names=True,filling_values =np.NaN) #fill with NaNs so must use their bounds
  Years=data["Year"][:]       ####<<<<<<<<<<<<<<<<<<<<<<<<
  
  RunName = "Tests17Aug_testing_E4"
  TraceFile = RunName + "_Trace.txt"
  tf=open(TraceFile,"w")
  print >>tf,"trace", RunName, datetime.datetime.now()
  tf.close()
  
  f = open(RunName+".txt","w")
  
  discoveryList={} 
  initialList = {}
  reps=1
  cutoff=0.5
  for i in range(reps):
    Rand = np.array([random.random() for y in Years])
#==============================================================================
#     Rand = np.array([ 0.81241285,  0.11963721,  0.81572623,  0.34903063,  0.56381217,
#         0.32776746,  0.56902457,  0.36930125,  0.18332472,  0.64712751,
#         0.29260813,  0.25061134,  0.94408559,  0.91931017,  0.72994937,
#         0.58337969,  0.63415929,  0.93789201,  0.19482533,  0.4010691 ,
#         0.28170177,  0.50161253,  0.42374253,  0.44809375,  0.27070208,
#         0.49046427,  0.19358991,  0.45453564,  0.47998151,  0.22703944,
#         0.73795753,  0.08065986,  0.14073268,  0.01014045,  0.3718777 ,
#         0.7655613 ,  0.78780508,  0.27818196,  0.06491864,  0.71508204,
#         0.55419559,  0.67285248,  0.24583931,  0.72373743,  0.78242748,
#         0.5162638 ,  0.60047823,  0.296487  ,  0.76244406,  0.46905657,
#         0.19566641,  0.23865389,  0.44041078,  0.10549116,  0.49024197,
#         0.60064108,  0.85073546,  0.00414606,  0.27871936,  0.82949543,
#         0.29738205,  0.76954792,  0.64555347,  0.83268187,  0.1355016 ,
#         0.28597666,  0.86491861,  0.53244471,  0.74712488,  0.06709093,
#         0.1068663 ,  0.09132201,  0.37564974,  0.38546621,  0.98983187,
#         0.97771759,  0.03422207,  0.95250432,  0.72937461,  0.9321408 ,
#         0.14136833,  0.26250987,  0.01194047,  0.11066821,  0.12297181,
#         0.6678355 ,  0.49281513,  0.74270458,  0.09798335,  0.06422087,
#         0.41264731,  0.1222206 ,  0.56102554,  0.11551947,  0.21535697,
#         0.01535822,  0.28325293,  0.403665  ,  0.9771023 ,  0.37017517,
#         0.61408612,  0.28172573,  0.11323753,  0.79806008,  0.83691362,
#         0.70245073,  0.59170967,  0.58994676,  0.10187166,  0.51172589,
#         0.21259711,  0.32266825,  0.07803303,  0.07952116,  0.68181968,
#         0.03778444,  0.8250724 ,  0.53384258,  0.282623  ,  0.20304002,
#         0.4744415 ,  0.05105567,  0.13277572,  0.4950562 ,  0.0692967 ,
#         0.97817368,  0.71669743,  0.72309605,  0.43852484,  0.83475738,
#         0.30544653,  0.81275605,  0.34389667,  0.50278543,  0.84740349,
#         0.051885  ,  0.76141463,  0.729362  ,  0.69606121,  0.10091101,
#         0.34650937,  0.9932162 ,  0.88042698,  0.67972356,  0.70877109,
#         0.49664535,  0.54708654,  0.73838778,  0.73180802,  0.65751174,
#         0.68796819,  0.74208573,  0.92626612,  0.62595156,  0.2055066 ,
#         0.721119  ,  0.58469861,  0.79470203,  0.37997864,  0.95067026,
#         0.81808836,  0.16648625,  0.24772792,  0.18402579,  0.6627768 ,
#         0.30631081,  0.67671502,  0.79034754,  0.03664392,  0.7263191 ,
#         0.31654487,  0.36966545,  0.41478369,  0.58828995,  0.21899307,
#         0.15285262,  0.41374726,  0.43519838,  0.90799403,  0.88727715,
#         0.03123549,  0.55225527,  0.68179498,  0.74810649,  0.97901864,
#         0.24103224,  0.27745015,  0.23412642,  0.53068976,  0.82463011,
#         0.70720961,  0.56664178,  0.99408574,  0.15750722,  0.01396173,
#         0.48261615,  0.41720631,  0.61062604,  0.21240388,  0.09021877,
#        0.37760615])
#==============================================================================
    '''
    Rand = np.array([ 0.64324326,  0.21930741,  0.1244528 ,  0.22169075,  0.84264838,
        0.67830664,  0.76266255,  0.25065429,  0.42688375,  0.97628721,
        0.12836764,  0.65160289,  0.60105608,  0.34408669,  0.04070487,
        0.40247888,  0.05195674,  0.92240849,  0.36342763,  0.43110751,
        0.65101805,  0.05058869,  0.22397402,  0.39195057,  0.18215362,
        0.23334128,  0.49453287,  0.86292751,  0.09639462,  0.15209089,
        0.33320689,  0.87971685,  0.87409072,  0.86387489,  0.96487674,
        0.17949379,  0.41827399,  0.63260275,  0.22113485,  0.04707216,
        0.55787911,  0.46454353,  0.57671954,  0.51742044,  0.48347653,
        0.90061547,  0.52892456,  0.67142054,  0.87360039,  0.61047586,
        0.85583316,  0.40558904,  0.36275935,  0.6536887 ,  0.02704021,
        0.48621901,  0.89989951,  0.68233809,  0.77150077,  0.17548277,
        0.43260673,  0.6641197 ,  0.9090032 ,  0.19756605,  0.4094763 ,
        0.62233147,  0.66959683,  0.25785318,  0.99114566,  0.77315398,
        0.33261257,  0.57552639,  0.59787583,  0.57309057,  0.28975642,
        0.78728485,  0.69851528,  0.67391941,  0.19099287,  0.30464899,
        0.65451396,  0.48135452,  0.01344684,  0.9755047 ,  0.01181421,
        0.98368475,  0.72510904,  0.91078699,  0.13711129,  0.44271504,
        0.06218034,  0.2946351 ,  0.4171683 ,  0.78694615,  0.80198248,
        0.19734225,  0.3197144 ,  0.29565569,  0.82268319,  0.59910277,
        0.12016079,  0.98395959,  0.21734559,  0.3651133 ,  0.67491085,
        0.87381084,  0.47512161,  0.45451124,  0.25736263,  0.73271135,
        0.80189047,  0.60671215,  0.38774219,  0.00691815,  0.54062939,
        0.14345425,  0.65805034,  0.11094342,  0.55590282,  0.86713355,
        0.87764494,  0.00870472,  0.24705249,  0.44403531,  0.47453249,
        0.63835235,  0.31642773,  0.336832  ,  0.4176651 ,  0.59378574,
        0.10963615,  0.38952407,  0.08568285,  0.3359799 ,  0.58528267,
        0.79713414,  0.00140362,  0.16690939,  0.75983224,  0.69796575,
        0.35040553,  0.17973782,  0.73377536,  0.02625453,  0.42325562,
        0.93603307,  0.47629316,  0.84271691,  0.70863257,  0.90216801,
        0.73399547,  0.52313723,  0.48540689,  0.69934059,  0.34094876,
        0.0659576 ,  0.14062765,  0.88490481,  0.00192505,  0.97681841,
        0.82281746,  0.53133003,  0.92359786,  0.89827533,  0.64627265,
        0.9306581 ,  0.06986758,  0.27236723,  0.42945854,  0.13347879,
        0.75387688,  0.39324943,  0.78518891,  0.61720469,  0.19222311,
        0.37791963,  0.17764423,  0.06739783,  0.39348858,  0.36374551,
        0.33006454,  0.95862939,  0.97268299,  0.43024978,  0.02893368,
        0.34206924,  0.61491774,  0.77305967,  0.07182273,  0.7623944 ,
        0.17261387,  0.93135717,  0.35318897,  0.90888097,  0.39728778,
        0.00745894,  0.32972421,  0.06115299,  0.06527679,  0.40686014,
        0.84256621])
    '''
    rand = np.array([ 0.49523659,  0.45579681,  0.8293389 ,  0.57066741,  0.36698356,
        0.50494943,  0.11468755,  0.83137379,  0.2877465 ,  0.16759134,
        0.49370929,  0.58479292,  0.70024498,  0.72533876,  0.47437117,
        0.50276379,  0.50002558,  0.09329194,  0.66580957,  0.39185105,
        0.74280447,  0.17516961,  0.48724713,  0.25587019,  0.0862142 ,
        0.33519199,  0.24765888,  0.75920359,  0.55374855,  0.71518447,
        0.0617846 ,  0.66101208,  0.7601255 ,  0.46609504,  0.95659172,
        0.56258313,  0.29310052,  0.32180116,  0.86647723,  0.64528137,
        0.22666548,  0.67697997,  0.42150835,  0.29937381,  0.22802639,
        0.77852367,  0.18228814,  0.80786725,  0.86995535,  0.33008637,
        0.30334666,  0.54278366,  0.23539573,  0.25427193,  0.81929239,
        0.69157558,  0.92476371,  0.77688199,  0.63705373,  0.45710281,
        0.76171356,  0.69462342,  0.56556555,  0.3179908 ,  0.64689819,
        0.28263754,  0.99285022,  0.30522544,  0.46313112,  0.45475422,
        0.55362966,  0.34515266,  0.35620343,  0.10047016,  0.14268542,
        0.32630085,  0.08152804,  0.43488547,  0.78727176,  0.06492124,
        0.99870039,  0.95638488,  0.88026366,  0.30936875,  0.25062926,
        0.19722009,  0.17978572,  0.40877558,  0.96807326,  0.10018074,
        0.98625507,  0.3720486 ,  0.41656155,  0.68943248,  0.49704996,
        0.58586196,  0.92905701,  0.80521425,  0.570477  ,  0.75240445,
        0.67265386,  0.49700543,  0.58088289,  0.57057158,  0.6034693 ,
        0.39237329,  0.20898083,  0.48904419,  0.36415124,  0.1513515 ,
        0.85179983,  0.02039117,  0.89953129,  0.17078232,  0.62698205,
        0.05707841,  0.60815797,  0.62962262,  0.02796569,  0.65308467,
        0.15302907,  0.69502492,  0.75135297,  0.09212866,  0.12363537,
        0.93408178,  0.45714621,  0.9217765 ,  0.45208139,  0.81271267,
        0.09441264,  0.03759134,  0.80298276,  0.77741424,  0.18537056,
        0.93615681,  0.35253899,  0.27275695,  0.09306062,  0.20971768,
        0.17224665,  0.93656314,  0.11602232,  0.10722419,  0.6577062 ,
        0.53148936,  0.49297903,  0.58619263,  0.42728485,  0.77298222,
        0.60734247,  0.91405847,  0.10471911,  0.59139347,  0.46264389,
        0.11961349,  0.74694292,  0.58033938,  0.78933342,  0.56811362,
        0.48039535,  0.50645863,  0.52453572,  0.27544549,  0.05845534,
        0.66517125,  0.71443783,  0.25851136,  0.65048435,  0.97670241,
        0.87644122,  0.7629569 ,  0.17786076,  0.04805553,  0.64412747,
        0.36614629,  0.48176892,  0.78053255,  0.77722434,  0.60837888,
        0.79204986,  0.85507658,  0.33369285,  0.14854147,  0.23811625,
        0.39242432,  0.79760548,  0.67498513,  0.43850435,  0.55105412,
        0.65654428,  0.50297278,  0.38049045,  0.1264982 ,  0.88987082,
        0.55213132,  0.54408984,  0.68642229,  0.25780613,  0.43582019,
        0.28787876])

    print >>f,datetime.datetime.now()
    for p in ["E4"]:#data.dtype.names:
      model = p
      if not model in discoveryList:
        discoveryList[model]={}
      if not model  in rogerslist:
        rogerslist[model] = []
      if not model in initialList:
        initialList[model]={}
      print >>f, model, RunName
      print "MODEL", model, i
      
    
      lo,hi=valid(data[model][:])   #####<<<<<<<<<<<<<<<<<<<<<<<<<<
      #if model != "Year": 
      if model != "Year":
        try:
          Dat=data[model][lo:hi]#-sigmoid.sig4(Years[lo:hi],data[model][lo:hi])
        except:
          print "Cannot estimate sigmoid"
          Dat=data[model][lo:hi]
        brks=convergentBreaks(Dat, Rand[lo:hi], Years[lo:hi], model, mode="control")
        print  >>f,"MODEL",model, "FOUND", brks
        br=tuple(brks[1]) #immutable so can use as key

        plotbreaksSummary(Years[lo:hi], data[model][lo:hi], [], title=model+" "+os.path.basename(fn), 
                      subtitle="BVave:"+str(br), breakyears=None, breakyearsAve=br, offset=0)
    
  #print discoveryList
  #print collate(discoveryList, cutoff=cutoff, reps=reps)
  f.close()
#    