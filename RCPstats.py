# -*- coding: utf-8 -*-
"""
Created on Thu Aug 06 16:41:26 2015

@author: s4493222
"""

SVNRevision="$Revision$"

#script to perform statistics for Roger (see ticket #24).

#1 we need to open all named csv files and create a simple data base from them.

import csv
import statbreaks
import os
import glob
import sys
import numpy as np

fn="C:\\Users\\s4493222\\Documents\\abrupt\\4Roger_RCP\\RCP26_05aug.x.txt.csv"
datapath="C:\\Users\\s4493222\\Documents\\abrupt\\c_test\\"
tdir="C:\\Users\s4493222\\Documents\\abrupt\c_test\\r10i1p1_CNRM-CM5_RCP8.5_icmip5_tas_Amon_ens_rcp45to85_0-360E_-90-90N_n_su_179\\r10i1p1_CNRM-CM5_RCP8.5_icmip5_tas_Amon_ens_rcp45to85_0-360E_-90-90N_n_su_179_0.trace"
def getannual(finalfn,datapath="C:\\Users\\s4493222\\Documents\\abrupt\\c_test\\"):
  #assume *.final.csv
  direc=os.path.splitext(os.path.splitext(finalfn)[0])[0]
  anfn=datapath+direc+"\\"+direc+"*_0.trace"
  anfn=glob.glob(anfn)
  if len(anfn) <1:
    #print "cannot find ",anfn,finalfn
    #sys.exit()
    return None,None,None,None
    
  else:
    anfn=anfn[0]
    #print os.path.exists(anfn)
    sb=statbreaks.brkrpt(anfn)
    return sb.breaks(), sb.years(), sb.ys(), sb.segments()
  
def statsfor(breaks, years, ys,segments) :
  #find 5 year averages centered on 2006 and 2095
  lowyr=list(years).index(2006)
  hiyr=list(years).index(2095)
  lowy=np.mean(ys[lowyr-2:lowyr+3])
  hiy=np.mean(ys[hiyr-2:hiyr+3])
  return lowy, hiy
  
 
with open(fn, 'rb') as csvfile:
  lines = csv.DictReader(csvfile, delimiter=',', quotechar='"')
  #print lines.fieldnames
  fields = ['Seq', 'Source', 'BreakSet', 'DatasetNo', 'SetNo', 'SetNoInData', 'Type', 'Model', 'Treat', '%BreakSet', 'Year', 'Predecessor', 'Successor', '%SegmentFound', 'Ti', 'Shift', '%YearOverall', 'Modal Year', 'ModalValue', 'SecondModalYear', 'SecondModalValue', 'SegmentTrend', 'PreTrend', 'PostTrend"', 'ShiftFromTrends', '30YearTrend', '15YrTi', '15YearShift', '30Yr Mode', '30YrModePct', '30YrMode2', '30YrMode2%', 'Stability', ':', 'ModelGroup', 'GroupCode', 'GroupNum', 'Model', 'GroupModNum', 'ModelNum', 'EnsembleCode', 'EnsembleNum', 'Scenario', ':', 'Byear', 'slope1', 'Pr(slope1)', 'Sig(slope1)', 'Rsq(slope1)', 'slope2', 'Pr(slope2)', 'Sig(slope2)', 'Rsq(slope2)', 'non-trend-shift', 'Pr(non-trend-shift)', 'Sig(non-trend-shift)', ':', 'PrTrends', 'sig(Trends)', 'Pr(shifts)', 'sig(shifts)', 'Pr(Regime)"', 'sig(Regime)', 'Modal']
  selfields=['Seq', 'Source', 'Type', 'Model', 'Treat', 'Ti', 'Shift', 'SegmentTrend', 'PreTrend', 'PostTrend"', 'ShiftFromTrends', 'Byear', 'slope1', 'Pr(slope1)', 'Sig(slope1)', 'Rsq(slope1)', 'slope2', 'Pr(slope2)', 'Sig(slope2)', 'Rsq(slope2)', 'non-trend-shift', 'Pr(non-trend-shift)', 'Sig(non-trend-shift)', 'PrTrends', 'sig(Trends)', 'Pr(shifts)', 'sig(shifts)', 'Pr(Regime)"', 'sig(Regime)', 'Modal'] 
#print a heading  
  for field in selfields: 
    print field,
  for field in ["LoMean","HiMean","DeltaTrend","n","CumShift"]:
    print field,
  print
#and all records
  state=0
  for line in lines:
    if line['Modal'] == "Y" and float((str(line['Byear']).strip())) >=2006:
      for field in selfields:
        print line[field],
    
      breaks, years, ys,segments=getannual(line['Source'])  
      if breaks !=None:
        los,his=statsfor(breaks, years, ys,segments)
        print los,his,
        bryr=float((str(line['Byear']).strip())) 
        bryr = min(bryr, 2095)
        nyear=bryr-yr1
        trdelta+=trend*nyear
        tshift += float((str(line['non-trend-shift']).strip()))
        try:
          trend=float((str(line['PostTrend"']).strip()))
        except:
          print "PostTrend", line['PostTrend"']
          
        yr1 = bryr
        print trdelta, nyear, tshift
      else:
        state=0
        print "End Rec"
    else:
      tshift=0.0
      trdelta=0.0
      try:
        trend=float((str(line['PostTrend"']).strip()))
      except:
        trend=0.0
        #print "PreTrend", line['PreTrend']
        pass
      #trend=float((str(line['PreTrend']).strip()))              
      yr1=2006