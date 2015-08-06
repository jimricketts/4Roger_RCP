# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 14:33:38 2015

@author: s4493222
"""

#shuffleRpt
#code to load a trace file and perform shuffle tests on it.

import numpy as np
import statbreaks
import shuffle
import scipy.stats as stats
import bivariate_multi as bivariate
import ConvergentBreaks as convergent_breaks
import random
import glob
import sys
import regress
import datetime
import os
SVNRevision="$Revision: 307 $"
raise "Deprecated! Use Shuffle_Report"
#def shuffle_stat(xs1, ys1, Years, pos, iterations=100):

class shuffler(object):
  def __init__(self, fn):
    self.__filename = fn
    self.__data=statbreaks.brkrpt(fn)
    self.__yearsegs, self.__ysegs, self.__xsegs = self.__data.segments()
    
  def data(self):    
    return self.__data
  def filename(self):
    return self.__filename

  def spanned(self, brknumber, span=6):
    loxs=self.__xsegs[brknumber-1]
    loys=self.__ysegs[brknumber-1]
    loyrs=self.__yearsegs[brknumber-1]
    hixs=self.__xsegs[brknumber]
    hiys=self.__ysegs[brknumber]
    hiyrs=self.__yearsegs[brknumber]
    if span!=0:
      xs=list(loxs[-span:])
      xs.extend(hixs[:span])
      ys=list(loys[-span:])
      ys.extend(hiys[:span])
      yrs=list(loyrs[-span:])
      yrs.extend(hiyrs[:span])
    else:
      xs=list(loxs[:])
      xs.extend(hixs[:])
      ys=list(loys[:])
      ys.extend(hiys[:])
      yrs=list(loyrs[:])
      yrs.extend(hiyrs[:])
    return xs, ys, yrs
    
  def shuffle(self, brknumber, offset=0, span=6):
    loxs=self.__xsegs[brknumber-1]
    loys=self.__ysegs[brknumber-1]
    loyrs=self.__yearsegs[brknumber-1]
    hixs=self.__xsegs[brknumber]
    hiys=self.__ysegs[brknumber]
    hiyrs=self.__yearsegs[brknumber]
    if span!=0:
      xs=list(loxs[-span:])
      xs.extend(hixs[:span])
      ys=list(loys[-span:])
      ys.extend(hiys[:span])
      yrs=list(loyrs[-span:])
      yrs.extend(hiyrs[:span])
      stat=shuffle.shuffle_stat(xs, ys, yrs, span+offset)
    else:
      xs=list(loxs[:])
      xs.extend(hixs[:])
      ys=list(loys[:])
      ys.extend(hiys[:])
      yrs=list(loyrs[:])
      yrs.extend(hiyrs[:])
    
      stat=shuffle.shuffle_stat(xs, ys, yrs, len(loxs) + offset)
    bins=np.bincount(stat[0])
    mode=bins.argmax()
    bv=bivariate.bivariate(ys, xs, anomalise=False, pr=0.01)
    return stat, bins[yrs[0]:], yrs, mode, stats.norm.fit(stat[0]),stats.norm.fit(stat[1]),stats.norm.fit(stat[2]), bv.maxTi(), bv.stepChange()


def exemplars(path):
  allfiles=glob.glob(path)
  uniques={}
  for f in allfiles:
    s = shuffler(f)
    bs=str(s.data().breaks())
    if not bs in uniques:
      uniques[bs]=[]
    uniques[bs].append(f)
  return uniques

def evalbreaks(shuf):
  breaklist = {}
  yearlist = {}      
  yearseglist={}
 # analysis=r[:-2]
  for bs in shuf.keys():
    with open(shuf[bs][0]+'.shuffle', 'w') as outf:
      b=shuffler(shuf[bs][0])
      for line in shuf[bs]:
        print >>outf, line
      #bs=str(shuf.data().breaks())
      print >>outf,bs
      if not bs in breaklist:
        breaklist[bs]=0
      breaklist[bs] += 1  
      bse=eval(bs)
      for bry in bse:
        if not bry in yearlist:
          yearlist[bry] = 0
        yearlist[bry] += 1
      for i in range(1,len(bs) -1):
        ysl=str(bse[i-1:i+2])
        if not ysl in yearseglist:
          yearseglist[ysl] =0
        yearseglist[ysl] += 1
        
      for i in range(1,len(b.data().breaks())-1):
        for offset in [-2,-1,0,1,2]:
          st=b.shuffle(i,offset=offset)
          print >>outf,i, b.data().breaks()[i], offset,st[3:],
          histog=[(st[2][j],  st[1][j]) for j in range(len(st[1]))]
          print >>outf, histog
        
      for i in range(1,len(b.data().breaks())-1):
        xs,ys,yrs=b.spanned(i)
        pcb=convergent_breaks.resample_break(np.array(ys), np.array(yrs), N=100,withmode=True)
        print >>outf,"C",pcb#, detrend=False)
        print >>outf,"S",b.shuffle(i, span=0)[-2]
      
    

if __name__ == "__main__":
#  ex=exemplars("r1i1p2_GISS-E2-H_RCP4.5_icmip5_tas_Amon_ens_rcp45to85_0-360E_-90-90N_n_su_052_/r1i1p2_GISS-E2-H_RCP4.5_icmip5_tas_Amon_ens_rcp45to85_0-360E_-90-90N_n_su_052*.trace")
#  evalbreaks(ex)
#  a=shuffler(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\4Roger_Nature_SVN_264\\HadCRUT.4.2.0.0.annual_ns_avg//HadCRUT.4.2.0.0.annual_ns_avg.txt_0.trace")
#  for i in [1,2,3]:
#    st=a.shuffle(i)
#  
#    print st[1:]
## a simple test of linear trend and shifts.
#  f=open ("log.log.txt", "w")
#  print >>f,"shfact, trfact,ystd,shift, shdate, bv.maxIndexTi(), bv.maxTi(), bv.stepChange(), int(bv.maxTi()>=bv.critical())"
#  for shfact in [0., 1., 2., 3., 4.]:
#    for trfact in [0., .5, 1., 1.5, 2., 2.5, 3., 3.5, 4.]:
#      for j in range(100):
#        ys = np.array([random.random() for i in range(100)])
#        xs = np.array([random.random() for i in range(100)])
#        ystd=stats.norm.fit(ys)[1]
#        trend = trfact* ystd/100.
#        shift = shfact*ystd/2.
#        shdate = int((0.5+random.random()) * 50)
#        yshift= np.array([ int(i >= shdate) * shift for i in range(len(ys))])
#        ytrend = np.array([ i * trend for i in range(len(ys))])
#        #print yshift
#        bv=bivariate.bivariate(ys + ytrend + yshift, xs, anomalise=False, pr=0.01)
#        print >>f,shfact, trfact,ystd,shift, shdate, bv.maxIndexTi(), bv.maxTi(), bv.stepChange(), int(bv.maxTi()>=bv.critical())
#  f.close()
#==============================================================================
#   b=shuffler("C:\\Users\\s4493222\\Documents\\abrupt\\4Roger_Nature_SVN_264\\r1i1p1_ACCESS1-3_RCP8.5_icmip5_tas_Amon_ens_rcp45to85_0-360E_-90-90N_n_su_156\\r1i1p1_ACCESS1-3_RCP8.5_icmip5_tas_Amon_ens_rcp45to85_0-360E_-90-90N_n_su_156.dat_98.trace")
#   for trf in glob.glob(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\4Roger_Nature_SVN_264\\HadCRUT.4.2.0.0.annual_ns_avg//HadCRUT.4.2.0.0.annual_ns_avg.txt_*.trace"):
#   for trf in glob.glob(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\4Roger_Nature_SVN_264\\HadCRUT.4.2.0.0.annual_30S_30N//HadCRUT.4.2.0.0.annual_30S_30N.txt_*.trace"):
#==============================================================================
#  with open("runlist_orig.txt","r") as r:
#    runs=r.readlines()
##  with open("runlist.txt","r") as r1:
##    runs_done=r1.readlines()
#  runs=["r1i1p2_GISS-E2-H_RCP4.5_icmip5_tas_Amon_ens_rcp45to85_0-360E_-90-90N_n_su_052/\n",
#          "adjusted.csv_meanmodel/\n",
#          "HadCRUT.4.3.0.0.annual_sh/\n",
#  
#        ]
  runs=glob.glob("*/")
  for r in runs:#["had4_krig_annual_v2_0_0\\n"]:#runs:
    print r
    if True:#if not r in runs_done:
      breaklist={}
      yearlist={}
      yearseglist={}
      analysis=r[:-1]
      print analysis
      finalname=analysis+'.final.csv'
      for trf in glob.glob("./"+analysis+"//*.trace"):
        #print trf
        b=shuffler(trf)
        #[1850.0, 1929, 1978, 1996, 2014.0]
        with open(b.filename()+'.shuffle', 'w') as outf:
          bs=str(b.data().breaks())
          print >>outf,bs
          if not bs in breaklist:
            breaklist[bs]=0
          breaklist[bs] += 1  
          for bry in b.data().breaks():
            if not bry in yearlist:
              yearlist[bry] = 0
            yearlist[bry] += 1
          bse=eval(bs)
          for i in range(1,len(bs) -1):
            ysl=str(bse[i-1:i+2])
            if not ysl in yearseglist:
              yearseglist[ysl] =0
            yearseglist[ysl] += 1
            
          for i in range(1,len(b.data().breaks())-1):
            for offset in [-2,-1,0,1,2]:
              st=b.shuffle(i,offset=offset)
              print >>outf,i, b.data().breaks()[i], offset,st[3:],
              histog=[(st[2][j],  st[1][j]) for j in range(len(st[1]))]
              print >>outf, histog
            
          for i in range(1,len(b.data().breaks())-1):
            xs,ys,yrs=b.spanned(i)
            pcb=convergent_breaks.resample_break(np.array(ys), np.array(yrs), N=100,withmode=True)
            print >>outf,"C",pcb#, detrend=False)
            print >>outf,"S",b.shuffle(i, span=0)[-2]
      
      
  #      with open(b.filename()+'.shuffle', 'r') as inf:
  #        lines =inf.readlines()
  #        for line in lines:
  #          print line[:-1]
  
      outf=open(finalname,"w")  
      print finalname
      for brl in breaklist:
        bre=eval(brl)
        for Y in range(len(bre[1:-1])):
          
          Year=bre[1:-1][Y]
          lowyr= list(b.data().years()).index(bre[Y])
          highyr=list(b.data().years()).index(bre[Y+2])+1
          midyr=list(b.data().years()).index(bre[Y+1])+1
          
          print >>outf, '"'+str(brl)+'"', breaklist[brl], Year, bre[Y], bre[Y+2], yearseglist[str(bre[Y:Y+3])],
          pcb=convergent_breaks.resample_break(b.data().ys()[lowyr:highyr],b.data().years()[lowyr:highyr], N=100,withmode=True) 
          print >>outf,pcb[1][0], pcb[2][0], yearlist[Year], pcb[3][0][0], pcb[3][0][1],
          if pcb[3][1] == None:
            print >>outf,0, 0,
          else:
            print >>outf,pcb[3][1][0], pcb[3][1][1],
          betaseg, alphaseg=regress.regress(b.data().ys()[lowyr:highyr],b.data().years()[lowyr:highyr])
          betalow, alphalow=regress.regress(b.data().ys()[lowyr:midyr],b.data().years()[lowyr:midyr])
          betahigh, alphahigh=regress.regress(b.data().ys()[midyr:highyr],b.data().years()[midyr:highyr])
          print >>outf,betaseg, betalow, betahigh,
          
                
          #now do a diagnostic on 15 year enclosed
          low15 = max(0, midyr-8)
          hi15=midyr+8
          beta15, alpha15=regress.regress(b.data().ys()[low15:hi15],b.data().years()[low15:hi15])
          print >>outf,(alphahigh+betahigh * bre[Y+1])-(alphalow+betalow * bre[Y+1]), beta15, 
          pcb15=convergent_breaks.resample_break(b.data().ys()[low15:hi15],b.data().years()[low15:hi15], N=100,withmode=True) 
          print >>outf,pcb15[1][0], pcb15[2][0], pcb15[3][0][0], pcb15[3][0][1],
          if pcb15[3][1] == None:
            print >>outf,0, 0
          else:
            print >>outf,pcb15[3][1][0], pcb15[3][1][1]
      print >>outf, "closed", datetime.datetime.now().strftime(" %X,%a,%d-%b-%Y")
      print >>outf, "Version numbers"
      print >>outf, "shuffleRpt",SVNRevision
      print >>outf, "shuffle", shuffle.SVNRevision
      print >>outf, "bivariate_multi as bivariate" ,bivariate.SVNRevision
      print >>outf, "ConvergentBreaks as convergent_breaks",shuffle.SVNRevision
      print >>outf, "regress",regress.SVNRevision
      print >>outf, "statbreaks",statbreaks.SVNRevision
      outf.close()
      
