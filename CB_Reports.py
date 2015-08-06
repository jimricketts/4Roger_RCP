# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 13:41:22 2014

@author: s4493222
"""

#CB_Reports.py - take convergent_breaks objects and formats a report as a csv file with shuffle_test added
#convergentBreaks returns   initialBreaks, newbreaks, statlist
#statlist is a list of [(ystats, tstats, shiftstats)]
#where each of the three elements is (mean, stdev)
#e.g ([1850.0, 1920.0, 1996.0, 2012.0], [1850.0, 1920.0, 1996.0, 2012.0], [((1920.0, 0.0), (42.682625, 0.76312649), (0.20108237061765552, 0.0020284819355781795)), ((1996.0, 0.0), (52.760872, 0.6639939), (0.55426463964501127, 0.004255255432305689))])

#shuffle cut retuens   stats.norm.fit(TiPosIndex), stats.norm.fit(TiList),stats.norm.fit(ShiftList), float(bins[mode])/iterations, mode, 

#grep Returning *.trace|awk 'BEGIN {FS="->"}; {print $2}' |awk 'BEGIN {FS = "\\[\\("}; {print $1}' |sort|uniq -c
#MyDocuments/abrupt/4Roger_Nature_SVN_264/had4_krig_annual_v2_0_0
#grep Returning *.trace|awk 'BEGIN {FS="->"}; {print $2}' |awk 'BEGIN {FS = "\\[\\("}; {print $1}' |sort|uniq -c|sed 's/\[//g' |sed 's/\]//g'|sed 's/,//g'
# grep Returning *.trace|awk 'BEGIN {FS="->"}; {print $2}' |awk 'BEGIN {FS = "\\[\\("}; {print $1}' |sort|uniq -c|sed 's/\[//g' |sed 's/\]//g'|sed 's/,//g'|awk '{for (i=1;i<=NF;i++) {a[$i] +=$0;}} END {for (y in a) {print y " " a[y];}}'
# for d in `ls -d */`;do pushd $d; echo $d;grep Returning *.trace|awk 'BEGIN {FS="->"}; {print $2}' |awk 'BEGIN {FS = "\\[\\("}; {print $1}' |sort|uniq -c|sed 's/\[//g' |sed 's/\]//g'|sed 's/,//g'|awk -v dr=$d '{for (i=2;i<=NF;i++) {a[$i] +=$0;}} END {for (y in a) {print dr " " y " " a[y];}}';popd ;done >tabulated.rpt
import numpy as np
from CMIP3 import CMIP3gw
from CMIP5 import CMIP5gw
from NOAAascii import NOAAascii
from ICMP import ICMP
from HADCRUT4 import HADCRUT4
from HADCRUT4 import CRU
TRENDS=False
if TRENDS:
  import ConvergentBreaksTrends as convergent_breaks
else:
  import ConvergentBreaks as convergent_breaks  
import shuffle
import matplotlib.pyplot as plt
import datetime
import os
import glob
#import csv
import errno


SVNRevision="$Revision: 360 $"
SCREENPR=0.01


def mkdir_p(path):
    """ 'mkdir -p' in Python """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
def applyShuffle(cb, xs, ys, Years):
  breaks=[0]
  for b in cb[2]:
    #print b
    breaks.extend([list(Years[:]).index(round(b[0][0])) +1])
  breaks.append(len(Years))
  
  #print "BREAKS",breaks
  results = []
  for b in range(len(breaks))[1:-1]:
    if xs[breaks[b-1]:breaks[b+1]] == []:
      print "oops", breaks, b
    results.append(shuffle.shuffle_cut(xs[breaks[b-1]:breaks[b+1]], ys[breaks[b-1]:breaks[b+1]], Years[breaks[b-1]:breaks[b+1]], breaks[b]-breaks[b-1]))
  return results

def report(filename,cb, xs, ys, Years, HeaderComments):
  import bivariate_multi as bivariate
  sc=applyShuffle(cb, xs, ys, Years)
  init='"%s"' % cb[0] #inital bisection result
  refined='"%s"' % cb[1] #after refinement list of years
  with open(filename,"w") as out:
    for line in HeaderComments:
      print >>out,'"%s"' % (str(line),)
    print >>out, init, datetime.datetime.now().strftime(" %X,%a,%d-%b-%Y")
    print >>out, refined    
    print >>out,"BreakDate,CritTi, BreakMean, BreakStDev, BreakTi0, BreakTi0StDev, Shift, ShiftStd, ShuffledBreakdate, ShuffledBreakStDev, ShuffledBreakTi0, ShuffledBreakTi0StDev, ShuffledShift, ShuffledShiftStd, ShuffledYearHitRatio, ShuffledModeYear"
    i = 0
    for entry in cb[2]:
      scentry=sc[i]
      i += 1
      pr = bivariate.critTi(0.01, int(cb[1][i+1]-cb[1][i-1]))
      print >>out,'%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d' % \
      (int(round(cb[1][i])),pr, entry[0][0],entry[0][1],entry[1][0],entry[1][1],entry[2][0],entry[2][1], scentry[0][0],scentry[0][1],scentry[1][0],scentry[1][1],scentry[2][0],scentry[2][1],scentry[3], int(scentry[4] -1))

def gatherNOAANames(path):
  #"C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\ZonalLandSeaASCII_NCDC\\
  import glob
  import os
  fnames = glob.glob(path+"\\*.asc")
  return fnames

def gatherCMIP5Names(path):
  import glob
  import os
  fnames = glob.glob(path+"\\*.GW")
  histfnames=glob.glob(path+"\\*historical*.GW")
  #build a dict or historical file names
  filelist=[]
  for h in histfnames:
    fn = os.path.basename(h)
    (var,sort, gcm, _, rep, _)=fn.split("_")
    for rcp in ["rcp26", "rcp45", "rcp60", "rcp85"]:
      candidates=glob.glob(path+"\\"+"%s_%s_%s_%s_%s_*.GW" % (var, sort, gcm, rcp,rep))
      for c in candidates:
        if h != c:
          #print rcp, len(candidates), os.path.basename(h), os.path.basename(c)
          filelist.append([h,c])
  return filelist

def gatherCMIP5ZonalNames(path):
  import glob
  import os
  fnames = glob.glob(path+"\\*.ZW")
  histfnames=glob.glob(path+"\\*historical*.ZW")
  #build a dict or historical file names
  filelist=[]
  for h in histfnames:
    fn = os.path.basename(h)
    (var,sort, gcm, _, rep, _, lo,hi)=fn.split("_")
    for rcp in ["rcp26", "rcp45", "rcp60", "rcp85"]:
      candidates=glob.glob(path+"\\"+"%s_%s_%s_%s_%s_*_%s_%s" % (var, sort, gcm, rcp,rep,lo,hi))
      for c in candidates:
        if h != c:
          #print rcp, len(candidates), os.path.basename(h), os.path.basename(c)
          filelist.append([h,c])
  return filelist

def gatherCMIP3Names(path):
  import glob
  fnames = glob.glob(path+"\\*_gw.txt")
  #build a dict or historical file names
  return fnames
  
def gatherNCDCZonalNames(path):
  import glob
  fnames = glob.glob(path+"\\*asc")
  #build a dict or historical file names
  return fnames
  
def gatherNCDCcsvNames(path):
  import glob
  fnames = glob.glob(path+"\\*4.csv")
  #build a dict or historical file names
  return fnames

def collated(path, intype, initial_header=0):
  import glob
  import os
  #import csvfile
  outf=open(path+"\\collated.csv","w",0) #open with forced flush because of huger line lengths
  tracenames=glob.glob(path+"\\*.trace")
  initialfn=None
  #now open the csv file f
  
  i = 0
  while i in range(len(tracenames)) and initialfn == None: #find the initial file
    tn = tracenames[i]
   # bn=os.path.basename(tn)
    csvfn=os.path.splitext(tn)[0]+".csv"
    csvfn=csvfn.replace("\\\\","\\") #fudge to overcome an error I made earlier
    if intype in []:
      n = 3
    else:
      n = 4
    if intype==None:
      raise Exception("Specify type of input")
    try:
      analysis=csvfile.CSVfile(csvfn, timevar="BreakDate",skip_header=n)
      initialfn=analysis.header()[0].replace("\\\\","\\").replace("\n","").replace('"','').replace("'","")
      c3=csvfile.CSVfile(initialfn,skip_header=initial_header)
      Years=c3.data4name("Year")
    except:
      pass
    i +=1
    
  #now a reporting pass  
  rept={}
  for tn in tracenames:
    #bn=os.path.basename(tn)
    csvfn=os.path.splitext(tn)[0]+".csv"
    csvfn=csvfn.replace("\\\\","\\") #fudge to overcome an error I made earlier
    analysis=csvfile.CSVfile(csvfn, timevar="BreakDate",skip_header=n)
    model=analysis.header()[1].replace("\n","").replace('"','').replace("'","")
#    if model == "BCCA2r1":
#      print model
    ys=c3.data4name(model)
    Yrs=c3.time4name(model)
    Tis=np.zeros(np.shape(Yrs),dtype=np.float32)
    shifts=np.zeros(np.shape(Yrs),dtype=np.float32)
    bdates=analysis["BreakDate"]
    pltbreaks=[Yrs[0]]
    for i in range(len(bdates)):
      for y in range(len(Yrs)):
        if Yrs[y] == bdates[i]:
          Tis[y] = analysis["BreakTi0"][i]
          shifts[y]=analysis["Shift"][i]
          pltbreaks.append(Yrs[y])
    rept[model]=(Yrs,ys,Tis,shifts)
    pltbreaks.append(Yrs[-1])
    graph(ys, Yrs, pltbreaks,model,csvfn)
  print >>outf, "Year",
  headerfields=c3.fields()
  headerfields.sort()
  for model in headerfields:
    if model != "Year":
      print >>outf,","+model+",",model+"_Ti0,",model+"_shift",
  print >>outf
  keys=rept.keys()
  keys.sort()
  for y in range(len(Years)):
    print >>outf,"%d" % (Years[y],),
    for k in keys:
      res=rept[k]
      ry=res[0].searchsorted(Years[y])
#      if k == "BCCA2r1":
#        print k
      if ry in range(len(res[0])) and Years[y] == res[0][ry]: #not the case if search goes out of bounds
        try:
          print >>outf,",%f,%f,%f" % (res[1][ry], res[2][ry], res[3][ry]),
        except:
          print >>outf,", , ,",
          pass
      else:
        print >>outf,", , ,",
    print >>outf

  headerfields=c3.fields()
  headerfields.sort()
  extent=["","","_Ti0","_shift"]
  
  tanoms={}
  for k in headerfields:
    if k !="Year":
      res=rept[k]
      lo=res[0].searchsorted(1961)
      hi=res[0].searchsorted(1990)+1
      tanoms[k] = np.mean(res[1][lo:hi])
    
  for j in [1,2,3]:
    print >>outf, "Year",
    for model in headerfields:
      if model != "Year":
        print >>outf,","+model+extent[j],
    print >>outf
    keys=rept.keys()
    keys.sort()
    for y in range(len(Years)):
      print >>outf,"%d" % (Years[y],),
      for k in keys:
        res=rept[k]
        ry=res[0].searchsorted(Years[y])
  #      if k == "BCCA2r1":
  #        print k
        if ry in range(len(res[0])) and Years[y] == res[0][ry]: #not the case if search goes out of bounds
          try:
            if j == 1:
              print >>outf,",%f" % (res[j][ry]-tanoms[k], ),
            else:
              print >>outf,",%f" % (res[j][ry], ),
          except:
            print >>outf,",",
            pass
        else:
          print >>outf,",",
      print >>outf
  outf.close()
  
def collatedCMIP5(path, intype):
  import glob
  import os
  #import csvfile
  outf=open(path+"\\collated_"+intype+"_.csv","w",0) #open with forced flush because of huger line lengths
  tracenames=glob.glob(path+"\\*_"+str(intype)+"_*.trace")

  rept={}
  headerfields=[]
  minyr=10000
  maxyr=0
  for tn in tracenames:
    #bn=os.path.basename(tn)
    csvfn=os.path.splitext(tn)[0]+".csv"
    csvfn=csvfn.replace("\\\\","\\") #fudge to overcome an error I made earlier
    analysis=csvfile.CSVfile(csvfn, timevar="BreakDate",skip_header=4)
    fn=analysis.header()[:2]
    fn[0]=fn[0].replace("\n","").replace('"','').replace("'","")
    fn[1]=fn[1].replace("\n","").replace('"','').replace("'","")
    
    model=analysis.header()[1].replace("\n","").replace('"','').replace("'","")
    model=os.path.basename(model).split("_")
    model=model[2]+"_"+model[3]+"_"+model[4]
    c3=CMIP5gw(fn)
    ys =  c3.Warming()
    Yrs= c3.Years()
    minyr=min(minyr,np.min(Yrs))
    maxyr=min(2100,max(maxyr,np.max(Yrs)))
#    if model == "BCCA2r1":
#      print model
    Tis=np.zeros(np.shape(Yrs),dtype=np.float32)
    shifts=np.zeros(np.shape(Yrs),dtype=np.float32)
    bdates=analysis["BreakDate"]
    pltbreaks=[Yrs[0]]
    for i in range(len(bdates)):
      for y in range(len(Yrs)):
        if Yrs[y] == bdates[i]:
          Tis[y] = analysis["BreakTi0"][i]
          shifts[y]=analysis["Shift"][i]
          pltbreaks.append(Yrs[y])
    rept[model]=(Yrs,ys,Tis,shifts)
    pltbreaks.append(Yrs[-1])
    graph(ys, Yrs, pltbreaks,model,csvfn)
    headerfields.append(model)
  Years=range(int(round(minyr)), int(round(maxyr)) + 1)
  print >>outf, "Year",

  headerfields.sort()
  for model in headerfields:
    if model != "Year":
      print >>outf,","+model+",",model+"_Ti0,",model+"_shift",
  print >>outf
  
  keys=headerfields
  for y in range(len(Years)):
    print >>outf,"%d" % (Years[y],),
    for k in keys:
      res=rept[k]
      ry=res[0].searchsorted(Years[y])
#      if k == "BCCA2r1":
#        print k
      if ry in range(len(res[0])) and Years[y] == res[0][ry]: #not the case if search goes out of bounds
        try:
          print >>outf,",%f,%f,%f" % (res[1][ry], res[2][ry], res[3][ry]),
        except:
          print >>outf,", , ,",
          pass
      else:
        print >>outf,", , ,",
    print >>outf

#  headerfields=c3.fields()
#  headerfields.sort()
  extent=["","","_Ti0","_shift"]
  
  tanoms={}
  for k in headerfields:
    if k !="Year":
      res=rept[k]
      lo=res[0].searchsorted(1961)
      hi=res[0].searchsorted(1990)+1
      tanoms[k] = np.mean(res[1][lo:hi])
    
  for j in [1,2,3]:
    print >>outf, "Year",
    for model in headerfields:
      if model != "Year":
        print >>outf,","+model+extent[j],
    print >>outf
    for y in range(len(Years)):
      print >>outf,"%d" % (Years[y],),
      for k in keys:
        res=rept[k]
        ry=res[0].searchsorted(Years[y])
  #      if k == "BCCA2r1":
  #        print k
        if ry in range(len(res[0])) and Years[y] == res[0][ry]: #not the case if search goes out of bounds
          try:
            if j == 1:
              print >>outf,",%f" % (res[j][ry]-tanoms[k], ),
            else:
              print >>outf,",%f" % (res[j][ry], ),
          except:
            print >>outf,",",
            pass
        else:
          print >>outf,",",
      print >>outf
  outf.close()

def collatedNCDC(path, intype, part):
  import glob
  import os
  #import csvfile
  outf=open(path+"\\"+part+"_collated.csv","w",0) #open with forced flush because of huger line lengths
  tracenames=glob.glob(path+"\\"+part+"*1880-2014.csv*.trace")
  initialfn=None
  #now open the csv file f
  
  i = 0
  while i in range(len(tracenames)) and initialfn == None: #find the initial file
    tn = tracenames[i]
   # bn=os.path.basename(tn)
    csvfn=os.path.splitext(tn)[0]+".csv"
    csvfn=csvfn.replace("\\\\","\\") #fudge to overcome an error I made earlier
    if intype in [0]:
      n = 3
    else:
      n = 4
    if intype==None:
      raise Exception("Specify type of input")
    try:
      analysis=csvfile.CSVfile(csvfn, timevar="BreakDate",skip_header=n)
      initialfn=analysis.header()[0].replace("\\\\","\\").replace("\n","").replace('"','').replace("'","")
      c3=csvfile.CSVfile(initialfn,skip_header=2)
      Years=c3.data4name("Year")
    except:
      pass
    i +=1
    
  #now a reporting pass  
  rept={}
  headerfields=[]
  for tn in tracenames:
    #bn=os.path.basename(tn)
    csvfn=os.path.splitext(tn)[0]+".csv"
    csvfn=csvfn.replace("\\\\","\\") #fudge to overcome an error I made earlier
    analysis=csvfile.CSVfile(csvfn, timevar="BreakDate",skip_header=n)
    initialfn=analysis.header()[0].replace("\\\\","\\").replace("\n","").replace('"','').replace("'","")
    c3=csvfile.CSVfile(initialfn,skip_header=2)
    model=os.path.basename(analysis.header()[0].replace("\n","").replace('"','').replace("'","")).split('.')[0]
    model=csvfn.split("_")[-1]
    headerfields.append(model)
#    if model == "BCCA2r1":
#      print model
    
    ys=c3.data4name("Value")
    Yrs=c3.time4name("Value")
    Tis=np.zeros(np.shape(Yrs),dtype=np.float32)
    shifts=np.zeros(np.shape(Yrs),dtype=np.float32)
    bdates=analysis["BreakDate"]
    pltbreaks=[Yrs[0]]
    for i in range(len(bdates)):
      for y in range(len(Yrs)):
        if Yrs[y] == bdates[i]:
          Tis[y] = analysis["BreakTi0"][i]
          shifts[y]=analysis["Shift"][i]
          pltbreaks.append(Yrs[y])
    rept[model]=(Yrs,ys,Tis,shifts)
    pltbreaks.append(Yrs[-1])
    graph(ys, Yrs, pltbreaks,model,csvfn)
  print >>outf, "Year",
  headerfields.sort()
  for model in headerfields:
    if model != "Year":
      print >>outf,","+model+",",model+"_Ti0,",model+"_shift",
  print >>outf
  keys=rept.keys()
  keys.sort()
  for y in range(len(Years)):
    print >>outf,"%d" % (Years[y],),
    for k in keys:
      res=rept[k]
      ry=res[0].searchsorted(Years[y])
#      if k == "BCCA2r1":
#        print k
      if ry in range(len(res[0])) and Years[y] == res[0][ry]: #not the case if search goes out of bounds
        try:
          print >>outf,",%f,%f,%f" % (res[1][ry], res[2][ry], res[3][ry]),
        except:
          print >>outf,", , ,",
          pass
      else:
        print >>outf,", , ,",
    print >>outf

  headerfields.sort()
  extent=["","","_Ti0","_shift"]
  
  tanoms={}
  for k in headerfields:
    if k !="Year":
      res=rept[k]
      lo=res[0].searchsorted(1961)
      hi=res[0].searchsorted(1990)+1
      tanoms[k] = np.mean(res[1][lo:hi])
    
  for j in [1,2,3]:
    print >>outf, "Year",
    for model in headerfields:
      if model != "Year":
        print >>outf,","+model+extent[j],
    print >>outf
    keys=rept.keys()
    keys.sort()
    for y in range(len(Years)):
      print >>outf,"%d" % (Years[y],),
      for k in keys:
        res=rept[k]
        ry=res[0].searchsorted(Years[y])
  #      if k == "BCCA2r1":
  #        print k
        if ry in range(len(res[0])) and Years[y] == res[0][ry]: #not the case if search goes out of bounds
          try:
            if j == 1:
              print >>outf,",%f" % (res[j][ry]-tanoms[k], ),
            else:
              print >>outf,",%f" % (res[j][ry], ),
          except:
            print >>outf,",",
            pass
        else:
          print >>outf,",",
      print >>outf
  outf.close()

def collatedNOAA(path, intype):
#  for fn in files:
#    print fn
#    data=NOAAascii(fn).annually()
#    #print data.annually()
#    ys=np.array([row[2] for row in data]) #2014 not complete
#    
#    Years=np.array([row[0] for row in data])
  import glob
  import os
  #import csvfile
  outf=open(path+"\\collatedNOAA.csv","w",0) #open with forced flush because of huger line lengths
  tracenames=glob.glob(path+"\\*asc.trace")
  initialfn=None
  #now open the csv file f
  
  #now a reporting pass  
  rept={}
  headerfields=[]
  for tn in tracenames:
    #bn=os.path.basename(tn)
    csvfn=os.path.splitext(tn)[0]+".csv"
    csvfn=csvfn.replace("\\\\","\\") #fudge to overcome an error I made earlier
    #print csvfn
    try:
      analysis=csvfile.CSVfile(csvfn, timevar="BreakDate",skip_header=3)
    except:
      analysis=csvfile.CSVfile(csvfn, timevar="BreakDate",skip_header=3)
    initialfn=analysis.header()[0].replace("\\\\","\\").replace("\n","").replace('"','').replace("'","")
    data=NOAAascii(initialfn).annually()
#    #print data.annually()
    ys=np.array([row[2] for row in data]) #2014 not complete
#    
    Yrs=np.array([row[0] for row in data]) 
    Years=Yrs
    model=os.path.basename(analysis.header()[0].replace("\n","").replace('"','').replace("'","")).split(".")
    model=(model[3]+"."+model[4]+"_"+model[5]+"."+model[6]).split("_")[0]
    print model
    headerfields.append(model)
#    if model == "BCCA2r1":
#      print model
    
    Tis=np.zeros(np.shape(Yrs),dtype=np.float32)
    shifts=np.zeros(np.shape(Yrs),dtype=np.float32)
    bdates=analysis["BreakDate"]
    pltbreaks=[Yrs[0]]
    for i in range(len(bdates)):
      for y in range(len(Yrs)):
        if Yrs[y] == bdates[i]:
          Tis[y] = analysis["BreakTi0"][i]
          shifts[y]=analysis["Shift"][i]
          pltbreaks.append(Yrs[y])
    rept[model]=(Yrs,ys,Tis,shifts)
    pltbreaks.append(Yrs[-1])
    graph(ys, Yrs, pltbreaks,model,csvfn)
  print >>outf, "Year",
  
  headerfields.sort()
  for model in headerfields:
    if model != "Year":
      print >>outf,","+model+",",model+"_Ti0,",model+"_shift",
  print >>outf
  keys=headerfields
  for y in range(len(Years)):
    print >>outf,"%d" % (Years[y],),
    for k in keys:
      res=rept[k]
      ry=res[0].searchsorted(Years[y])
#      if k == "BCCA2r1":
#        print k
      if ry in range(len(res[0])) and Years[y] == res[0][ry]: #not the case if search goes out of bounds
        try:
          print >>outf,",%f,%f,%f" % (res[1][ry], res[2][ry], res[3][ry]),
        except:
          print >>outf,", , ,",
          pass
      else:
        print >>outf,", , ,",
    print >>outf

  extent=["","","_Ti0","_shift"]
  
  tanoms={}
  for k in headerfields:
    if k !="Year":
      res=rept[k]
      lo=res[0].searchsorted(1961)
      hi=res[0].searchsorted(1990)+1
      tanoms[k] = np.mean(res[1][lo:hi])
    
  for j in [1,2,3]:
    print >>outf, "Year",
    for model in headerfields:
      if model != "Year":
        print >>outf,","+model+extent[j],
    print >>outf
    keys=rept.keys()
    keys.sort()
    for y in range(len(Years)):
      print >>outf,"%d" % (Years[y],),
      for k in keys:
        res=rept[k]
        ry=res[0].searchsorted(Years[y])
  #      if k == "BCCA2r1":
  #        print k
        if ry in range(len(res[0])) and Years[y] == res[0][ry]: #not the case if search goes out of bounds
          try:
            if j == 1:
              print >>outf,",%f" % (res[j][ry]-tanoms[k], ),
            else:
              print >>outf,",%f" % (res[j][ry], ),
          except:
            print >>outf,",",
            pass
        else:
          print >>outf,",",
      print >>outf
  outf.close()

def collatedICMP(path, intype, initial_header=0):
  import glob
  import os
  #import csvfile
  outf=open(path+"\\collated.csv","w",0) #open with forced flush because of huger line lengths
  tracenames=glob.glob(path+"\\*.trace")
  initialfn=None
  #now open the csv file f
  
  i = 0
  while i in range(len(tracenames)) and initialfn == None: #find the initial file
    tn = tracenames[i]
   # bn=os.path.basename(tn)
    csvfn=os.path.splitext(tn)[0]+".csv"
    csvfn=csvfn.replace("\\\\","\\") #fudge to overcome an error I made earlier
    if intype in []:
      n = 3
    else:
      n = 4
    if intype==None:
      raise Exception("Specify type of input")
    try:
      analysis=csvfile.CSVfile(csvfn, timevar="BreakDate",skip_header=n)
      initialfn=analysis.header()[0].replace("\\\\","\\").replace("\n","").replace('"','').replace("'","")
      c3=csvfile.CSVfile(initialfn,skip_header=initial_header)
      Years=c3.data4name("Year")
    except:
      pass
    i +=1
    
  #now a reporting pass  
  rept={}
  for tn in tracenames:
    #bn=os.path.basename(tn)
    csvfn=os.path.splitext(tn)[0]+".csv"
    csvfn=csvfn.replace("\\\\","\\") #fudge to overcome an error I made earlier
    analysis=csvfile.CSVfile(csvfn, timevar="BreakDate",skip_header=n)
    model=analysis.header()[1].replace("\n","").replace('"','').replace("'","")
#    if model == "BCCA2r1":
#      print model
    ys=c3.data4name(model)
    Yrs=c3.time4name(model)
    Tis=np.zeros(np.shape(Yrs),dtype=np.float32)
    shifts=np.zeros(np.shape(Yrs),dtype=np.float32)
    bdates=analysis["BreakDate"]
    pltbreaks=[Yrs[0]]
    for i in range(len(bdates)):
      for y in range(len(Yrs)):
        if Yrs[y] == bdates[i]:
          Tis[y] = analysis["BreakTi0"][i]
          shifts[y]=analysis["Shift"][i]
          pltbreaks.append(Yrs[y])
    rept[model]=(Yrs,ys,Tis,shifts)
    pltbreaks.append(Yrs[-1])
    graph(ys, Yrs, pltbreaks,model,csvfn)
  print >>outf, "Year",
  headerfields=c3.fields()
  headerfields.sort()
  for model in headerfields:
    if model != "Year":
      print >>outf,","+model+",",model+"_Ti0,",model+"_shift",
  print >>outf
  keys=rept.keys()
  keys.sort()
  for y in range(len(Years)):
    print >>outf,"%d" % (Years[y],),
    for k in keys:
      res=rept[k]
      ry=res[0].searchsorted(Years[y])
#      if k == "BCCA2r1":
#        print k
      if ry in range(len(res[0])) and Years[y] == res[0][ry]: #not the case if search goes out of bounds
        try:
          print >>outf,",%f,%f,%f" % (res[1][ry], res[2][ry], res[3][ry]),
        except:
          print >>outf,", , ,",
          pass
      else:
        print >>outf,", , ,",
    print >>outf

  headerfields=c3.fields()
  headerfields.sort()
  extent=["","","_Ti0","_shift"]
  
  tanoms={}
  for k in headerfields:
    if k !="Year":
      res=rept[k]
      lo=res[0].searchsorted(1961)
      hi=res[0].searchsorted(1990)+1
      tanoms[k] = np.mean(res[1][lo:hi])
    
  for j in [1,2,3]:
    print >>outf, "Year",
    for model in headerfields:
      if model != "Year":
        print >>outf,","+model+extent[j],
    print >>outf
    keys=rept.keys()
    keys.sort()
    for y in range(len(Years)):
      print >>outf,"%d" % (Years[y],),
      for k in keys:
        res=rept[k]
        ry=res[0].searchsorted(Years[y])
  #      if k == "BCCA2r1":
  #        print k
        if ry in range(len(res[0])) and Years[y] == res[0][ry]: #not the case if search goes out of bounds
          try:
            if j == 1:
              print >>outf,",%f" % (res[j][ry]-tanoms[k], ),
            else:
              print >>outf,",%f" % (res[j][ry], ),
          except:
            print >>outf,",",
            pass
        else:
          print >>outf,",",
      print >>outf
  outf.close()

def graph(ys, Years, breaks,title, savename):
  import regress
  fig=plt.figure()
  plt.plot(Years, ys)

  segments=[Years.searchsorted(b+1) for b in breaks]
  #The above is based on getting the break into the start of the next period - but this leaves the first point out by one
  segments[0]-=1
  lohi=[(segments[i], segments[i+1]+1) for i in range(len(segments)-1)]
  for lh in lohi:
    (lo, hi) = lh
    stats=regress.analysed_regress(ys[lo:hi], Years[lo:hi])
    yhat, _=regress.residuals(ys[lo:hi], Years[lo:hi], stats)
    plt.plot(Years[lo:hi], yhat,'r-')
    ym=[np.mean(yhat) for y in yhat]
    plt.plot(Years[lo:hi], ym, 'k-')
  fig.suptitle(title+"\n"+datetime.datetime.now().strftime("%X,%a,%d-%b-%Y "))
  #ax.set_title(title)
  #plt.show()
  if savename != None:
    fig.savefig(savename+'.png')
  plt.close(fig)
  

if __name__ == "__main__":
  import random
  import csvfile
  import UAH_MSU
  import RSS
  import sys
#UAH MSU data
#  for fn in UAH_MSU.Headers[2:]:# ["Globe"]:
#    for layer in ['tlt', 'tmt','ttp', 'tls']:
#      uah=UAH_MSU.UAH_MSU(layer)
#      dirname = os.path.basename(uah.filename())+"_"+fn 
#      newdir = dirname
#      mkdir_p(newdir)
##       fn1=os.path.basename(c3.filename())
#      picfilename=dirname+"_pic.png"
#      if os.path.exists(picfilename):
#        print picfilename,"exists!"
#      else:
#        ys=uah.annual(fn)
#        Years=np.array(range(1979,2015))  
#        counts={}
#        for iteration in range(100):
#          xs = np.array([random.random() for y in ys])
#          fn1=dirname+"//"+dirname+"_"+str(iteration)
#          convergent_breaks.TraceFile = fn1+'.trace'
#          cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#          report(fn1+'.csv',cb, xs, ys, Years, [fn1,fn])
#          breaks=cb[1]
#          sbreaks=str(breaks)
#          if not sbreaks in counts:
#            counts[sbreaks] = 0
#          counts[sbreaks] += 1
#        maxc=0
#        maxb=''
#        for b in counts.keys():
#          if maxb=='':
#            maxb=b
#            maxc=counts[b]
#          else:
#            if counts[b] > maxc:
#              maxb=b
#              maxc=counts[b]
#        graph(ys, Years, eval(maxb),fn1+" "+str(maxc), dirname+"_pic")
#        plt.close("all")
##RSS MSU data
##RSSMSUCovers
##  for fn in RSS.RSSBands:# ["Globe"]:
##   for layer in RSS.RSSdataset:
#  for cover in RSS.RSSMSUCovers:
#    for layer in RSS.RSSdataset:
#      rss=RSS.RSS_MSU(layer,cover)
#      for fn in rss.headers():# ["Globe"]:
#        dirname = os.path.basename(rss.filename(fn))+"_"+layer 
#        newdir = dirname
#        mkdir_p(newdir)
#  #       fn1=os.path.basename(c3.filename())
#        picfilename=dirname+"_pic.png"
#        if os.path.exists(picfilename):
#          print picfilename,"exists!"
#        else:
#          ys=rss.annual(fn)
#          Years=np.array(rss.yearrange() )
#          counts={}
#          for iteration in range(100):
#            xs = np.array([random.random() for y in ys])
#            fn1=dirname+"//"+dirname+"_"+str(iteration)
#            convergent_breaks.TraceFile = fn1+'.trace'
#            cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#            report(fn1+'.csv',cb, xs, ys, Years, [fn1,fn])
#            breaks=cb[1]
#            sbreaks=str(breaks)
#            if not sbreaks in counts:
#              counts[sbreaks] = 0
#            counts[sbreaks] += 1
#          maxc=0
#          maxb=''
#          for b in counts.keys():
#            if maxb=='':
#              maxb=b
#              maxc=counts[b]
#            else:
#              if counts[b] > maxc:
#                maxb=b
#                maxc=counts[b]
#          graph(ys, Years, eval(maxb),fn1+" "+str(maxc), dirname+"_pic")
#          plt.close("all")
#  sys.exit()
#@another test case' Hadley and CRU HADCRUT4 are subtly different
#==============================================================================
#  c3=csvfile.CSVfile("C:\\Users\\s4493222\\Documents\\abrupt\\Hadly_Vs_CRU_Versions_Of_HADCRUT4_hemispheric\\DifferencePlot.csv")
#  c3=csvfile.CSVfile("C:\\Users\\s4493222\\Documents\\abrupt\\Hadly_Vs_CRU_Versions_Of_HADCRUT4_hemispheric\\DifferencePlotSST.csv")
#  for fn in c3.fields():
#    if fn != "Year":
#      dirname = os.path.basename(c3.filename())+"_"+fn 
#      newdir = dirname
#      mkdir_p(newdir)
##       fn1=os.path.basename(c3.filename())
#      picfilename=dirname+"_pic.png"
#      if os.path.exists(picfilename):
#        print picfilename,"exists!"
#      else:
#        ys=c3.data4name(fn)
#        Years=c3.time4name(fn)       
#        counts={}
#        for iteration in range(100):
#          xs = np.array([random.random() for y in ys])
#          fn1=dirname+"//"+dirname+"_"+str(iteration)
#          convergent_breaks.TraceFile = fn1+'.trace'
#          cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#          report(fn1+'.csv',cb, xs, ys, Years, [fn1,fn])
#          breaks=cb[1]
#          sbreaks=str(breaks)
#          if not sbreaks in counts:
#            counts[sbreaks] = 0
#          counts[sbreaks] += 1
#        maxc=0
#        maxb=''
#        for b in counts.keys():
#          if maxb=='':
#            maxb=b
#            maxc=counts[b]
#          else:
#            if counts[b] > maxc:
#              maxb=b
#              maxc=counts[b]
#        graph(ys, Years, eval(maxb),fn1+" "+str(maxc), dirname+"_pic")
#        plt.close("all")
  
#@a test case
#  fn1=[
#  "C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_CSIRO-Mk3-6-0_historical_r10i1p1_185001-200512.GW",
#  "C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_CSIRO-Mk3-6-0_rcp45_r10i1p1_200601-210012.GW"
#  ]  
#  counts={}
#  newdir= "_tests"
#  mkdir_p(newdir)
#  icmp=CMIP5gw(fn1)
#  ys =  icmp.Warming()
#  xs = np.array([random.random() for y in ys])
#  Years= icmp.Years()
#  fn = newdir+"\\"+os.path.basename(fn1[1])
#  for iteration in range(100):
#    xs = np.array([random.random() for y in ys])
#    convergent_breaks.TraceFile = fn+"_"+str(iteration)+'.trace'
#    cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#    report(fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [fn1[1],fn])
#    breaks=cb[1]
#    sbreaks=str(breaks)
#    if not sbreaks in counts:
#      counts[sbreaks] = 0
#    counts[sbreaks] += 1
#  maxc=0
#  maxb=''
#  for b in counts.keys():
#    if maxb=='':
#      maxb=b
#      maxc=counts[b]
#    else:
#      if counts[b] > maxc:
#        maxb=b
#        maxc=counts[b]
#  graph(ys, Years, eval(maxb),os.path.basename(fn)+" "+str(maxc), os.path.basename(fn)+"_pic")
#  plt.close("all")
#  raise Exception("done")
#  c3=csvfile.CSVfile(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\CMIP5_breakpoints\\artifishul\\Artyfishul_dartar4Jim (1).csv")
#  for fn in c3.fields():
#    if fn != "Year":
#      ys=c3.data4name(fn)
#      Years=c3.time4name(fn)
#      xs = np.array([random.random() for y in ys])
#      convergent_breaks.TraceFile = c3.filename()+"_"+fn+'.trace'
#      cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#      report(c3.filename()+"_"+fn+'.csv',cb, xs, ys, Years, [c3.filename(),fn])
#
#  collated(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\CMIP5_breakpoints\\artifishul\\",0,initial_header=4)      
  
#  fn=[
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_MIROC5_historical_r1i1p1_185001-201212.GW',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_MIROC5_rcp85_r1i1p1_200601-210012.GW'
#  ]
#  bfn=[
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_NorESM1-M_historical_r1i1p1_185001-200512.GW',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_NorESM1-M_rcp85_r1i1p1_200601-210012.GW'
#  ]
#  
#==============================================================================
#   c3=csvfile.CSVfile("C:\\Users\\s4493222\\Documents\\ReferenceData\\Rogers_Analysis_17Jul2014\\\Data 4 Jim\\Global temps.csv")
#   for fn in c3.fields():
#     if fn != "Year":
#       ys=c3.data4name(fn)
#       Years=c3.time4name(fn)
#       xs = np.array([random.random() for y in ys])
#       convergent_breaks.TraceFile = c3.filename()+"_"+fn+'.trace'
#       cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#       report(c3.filename()+"_"+fn+'.csv',cb, xs, ys, Years, [c3.filename(),fn])
# 
#   collated("C:\\Users\\s4493222\\Documents\\ReferenceData\\Rogers_Analysis_17Jul2014\\Data 4 Jim",0)      
#==============================================================================




#  files=gatherCMIP5ZonalNames('C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\Zonal85')
#  for fn in files:
#  #for fn in [['C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_bcc-csm1-1-m_historical_r1i1p1_185001-201212.GW', 'C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_bcc-csm1-1-m_rcp26_r1i1p1_200601-210012.GW']]:
#  #for fn in [['C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_CSIRO-Mk3-6-0_historical_r4i1p1_185001-200512.GW', 'C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_CSIRO-Mk3-6-0_rcp26_r4i1p1_200601-210012.GW']]:
#    c3=CMIP5gw(fn)
#    ys =  c3.Warming()
#    xs = np.array([random.random() for y in ys])
#    Years= c3.Years()
#    convergent_breaks.TraceFile = fn[1]+'.trace'
#    cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#    report(fn[1]+'.csv',cb, xs, ys, Years, fn)
#  for intype in ["rcp85"]:
#    collatedCMIP5("C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\Zonal85",intype)      


#  files=gatherCMIP5Names('C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW')
#  for fn in files:
#  #for fn in [['C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_bcc-csm1-1-m_historical_r1i1p1_185001-201212.GW', 'C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_bcc-csm1-1-m_rcp26_r1i1p1_200601-210012.GW']]:
#  #for fn in [['C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_CSIRO-Mk3-6-0_historical_r4i1p1_185001-200512.GW', 'C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_CSIRO-Mk3-6-0_rcp26_r4i1p1_200601-210012.GW']]:
#    c3=CMIP5gw(fn)
#    ys =  c3.Warming()
#    xs = np.array([random.random() for y in ys])
#    Years= c3.Years()
#    convergent_breaks.TraceFile = fn[1]+'.trace'
#    cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#    report(fn[1]+'.csv',cb, xs, ys, Years, fn)
#  for intype in ["rcp26","rcp45","rcp60","rcp85"]:
#    collatedCMIP5("C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW",intype)      
#    
#==============================================================================
#   files =gatherCMIP3Names("C:\\Users\\s4493222\\Documents\\ReferenceData\\CSIRO_Model_Set\\gwFiles")
#   badf="C:\\Users\\s4493222\\Documents\\ReferenceData\\CSIRO_Model_Set\\gwFiles\\cnrm_cm3.sresa2.run1.monthly.tas_A1_1860-2099_gw.txt"
#   for fn in files:#[badf]:#files:
#     c3=CMIP3gw(fn)
#     ys =  c3.Warming()
#     xs = np.array([random.random() for y in ys])
#     Years= c3.Years()
#     try:
#       convergent_breaks.TraceFile = fn+'.trace'
#       cb=convergent_breaks.convergentBreaks(ys, xs, Years, os.path.basename(fn), mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#       
#     except Exception as e:
#       #convergent_breaks.tf.close()
#       print  os.path.basename(fn)
#       raise
#     report(fn+'.csv',cb, xs, ys, Years, [fn])
#==============================================================================
#==============================================================================
# 
#  import Global_1880_2014_18 as data
#  Years = data.data[:,3][:]
#  ys = data.data[:, 1][:]
#  xs = data.data[:, 4][:]
#  convergent_breaks.TraceFile = 'Global_1880_2014_18.trace'
#  cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#  report("rn_Global_1880_2014_18"+'.csv',cb, xs, ys, Years, ["Dummy"])
#==============================================================================
  
#==============================================================================

  ##files =gatherNCDCcsvNames("./Global_1880-2014.csv")
#==============================================================================
#  label = "SH"
#  for fn in ["./"+label+"_1880-2014.csv"]:
#    data=np.genfromtxt(fn ,skip_header=2, names=True, delimiter=",")
#    ys =  data["Value"]
#    Years= data["Year"]
#    for iteration in range(100):
#      if os.path.exists(fn+"__"+str(iteration)+'.trace'):
#        print iteration, fn+"__"+str(iteration)+'.trace exists!'
#      else:
#        xs = np.array([random.random() for y in ys])
#        convergent_breaks.TraceFile = fn+"__"+str(iteration)+'.trace'
#        cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#        report(fn+"__"+str(iteration)+'.csv',cb, xs, ys, Years, [fn])
#  collatedNCDC("./",0, label)
#==============================================================================

# 
#  files = gatherNOAANames("C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\ZonalLandSeaASCII_NCDC")
#  #fn="C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\ZonalLandSeaASCII_NCDC\\aravg.mon.land_ocean.90S.90N.v3.5.4.201408.asc"
#  for fn in files:
#    print fn
#    data=NOAAascii(fn).annually()
#    #print data.annually()
#    ys=np.array([row[2] for row in data]) #2014 not complete
#    
#    Years=np.array([row[0] for row in data])
#    xs = np.array([random.random() for y in ys])
#    convergent_breaks.TraceFile = fn+'.trace'
#    cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#    report(fn+'.csv',cb, xs, ys, Years, [fn])
#  collatedNOAA("C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\ZonalLandSeaASCII_NCDC",0) 
#  
##  c3=csvfile.CSVfile("C:\\Users\\s4493222\\Documents\\ReferenceData\\Rogers_Analysis_17Jul2014\\Data 4 Jim\\GISSTEMPv3\\GISSTEMP v3.csv",skip_header=8)
##  c3=csvfile.CSVfile("C:\\Users\\s4493222\\Documents\\ReferenceData\\Rogers_Analysis_17Jul2014\\Data 4 Jim\\GISSTEMPv3_2013\\GISSTEMPv3_2013.csv",skip_header=8)
#  c3=csvfile.CSVfile("GISSTEMPto6-2013b.csv",skip_header=9)
#  for fn in c3.fields():
#    if fn != "Year":
#      ys=c3.data4name(fn)
#      Years=c3.time4name(fn)
#      for iteration in range(100):
#        xs = np.array([random.random() for y in ys])
#        convergent_breaks.TraceFile = c3.filename()+"_"+fn+"_"+str(iteration)+'.trace'
#        cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#        report(c3.filename()+"_"+fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [c3.filename(),fn])
#
#  collated(".\\",0,initial_header=9)      

#######################################################################################################
#The following code for cmip 3 DATA FROM Roger in EXCEL files

#==============================================================================
#  c3=csvfile.CSVfile("C:\\Users\\s4493222\\Documents\\ReferenceData\\Rogers_Analysis_17Jul2014\\\Data 4 Jim\\Global temps.csv")
#  for fn in c3.fields():
#    if fn != "Year":
#      dirname = os.path.basename(c3.filename())+"_"+fn 
#      newdir = dirname
#      mkdir_p(newdir)
##       fn1=os.path.basename(c3.filename())
#      picfilename=dirname+"_pic.png"
#      if os.path.exists(picfilename):
#        print picfilename,"exists!"
#      else:
#        print fn
#        ys=c3.data4name(fn)
#        Years=c3.time4name(fn)       
#        counts={}
#        for iteration in range(100):
#          xs = np.array([random.random() for y in ys])
#          fn1=dirname+"//"+dirname+"_"+str(iteration)
#          convergent_breaks.TraceFile = fn1+'.trace'
#          cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#          report(fn1+'.csv',cb, xs, ys, Years, [fn1,fn])
#          breaks=cb[1]
#          sbreaks=str(breaks)
#          if not sbreaks in counts:
#            counts[sbreaks] = 0
#          counts[sbreaks] += 1
#        maxc=0
#        maxb=''
#        for b in counts.keys():
#          if maxb=='':
#            maxb=b
#            maxc=counts[b]
#          else:
#            if counts[b] > maxc:
#              maxb=b
#              maxc=counts[b]
#        graph(ys, Years, eval(maxb),fn1+" "+str(maxc), dirname+"_pic")
#        plt.close("all")
##==============================================================================
#
###==============================================================================
###NCDC Land, Ocean and Land/Ocean Hemispheric
###==============================================================================
#  from NCDCset import NCDCFileSet
#  fns=[
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_LandOcean\\SH_LandOcean*1880-2015.csv',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_LandOcean\\NH_LandOcean*1880-2015.csv',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_LandOcean\\Global_LandOcean*1880-2015.csv',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Land\\SH-Land*1880-2015.csv',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Land\\NH-Land*1880-2015.csv',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Land\\Global-Land*1880-2015.csv',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Ocean\\SH-Ocean*1880-2015.csv',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Ocean\\NH-Ocean*1880-2015.csv',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Ocean\\Global-Ocean*1880-2015.csv'
#  ]
#  for fn in fns:
#    data = NCDCFileSet(fn)
#    ys =  data.annual()
#    Years= data.years()
#    counts={}
#    newdir= os.path.splitext((os.path.basename(data.filename())))[0]
#    mkdir_p(newdir)
#    fn1=newdir+'//'+os.path.basename(data.filename())
#    picfn=os.path.basename(data.filename())+"_pic"  
#    if os.path.exists(picfn+".png"):
#      print picfn+".png", "exists!!"
#    else:
#      for iteration in range(100):
#        data = NCDCFileSet(fn)
#        xs = np.array([random.random() for y in ys])
#        convergent_breaks.TraceFile = fn1+"_"+str(iteration)+'.trace'
#        cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#        report(fn1+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [fn1])
#        breaks=cb[1]
#        sbreaks=str(breaks)
#        if not sbreaks in counts:
#          counts[sbreaks] = 0
#        counts[sbreaks] += 1
#      maxc=0
#      maxb=''
#      for b in counts.keys():
#        if maxb=='':
#          maxb=b
#          maxc=counts[b]
#        else:
#          if counts[b] > maxc:
#            maxb=b
#            maxc=counts[b]
#      graph(ys, Years, eval(maxb),fn1+"_"+str(iteration)+'.csv'+" "+str(maxc), picfn)
#      plt.close("all")
##
###==============================================================================
###NCDC Land, Ocean and Land/Ocean Hemispheric by month
###==============================================================================
#  from NCDCset import NCDCFileSet
#  fns=[
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_LandOcean\\SH_LandOcean*1880-2015.csv',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_LandOcean\\NH_LandOcean*1880-2015.csv',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_LandOcean\\Global_LandOcean*1880-2015.csv',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Land\\SH-Land*1880-2015.csv',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Land\\NH-Land*1880-2015.csv',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Land\\Global-Land*1880-2015.csv',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Ocean\\SH-Ocean*1880-2015.csv',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Ocean\\NH-Ocean*1880-2015.csv',
#  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Ocean\\Global-Ocean*1880-2015.csv'
#  ]
#  for fn in fns:
#    for month in range(12):    
#      data = NCDCFileSet(fn)
#      ys =  data[month]
#      Years= data.years()
#      counts={}
#      newdir= os.path.splitext((os.path.basename(data.filename(month))))[0]
#      mkdir_p(newdir)
#      fn1=newdir+'//'+os.path.basename(data.filename(month))
#      picfn=os.path.basename(data.filename(month))+"_pic"  
#      if os.path.exists(picfn+".png"):
#        print picfn+".png", "exists!!"
#      else:
#        for iteration in range(100):
#          data = NCDCFileSet(fn)
#          xs = np.array([random.random() for y in ys])
#          convergent_breaks.TraceFile = fn1+"_"+str(iteration)+'.trace'
#          cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#          report(fn1+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [fn1])
#          breaks=cb[1]
#          sbreaks=str(breaks)
#          if not sbreaks in counts:
#            counts[sbreaks] = 0
#          counts[sbreaks] += 1
#        maxc=0
#        maxb=''
#        for b in counts.keys():
#          if maxb=='':
#            maxb=b
#            maxc=counts[b]
#          else:
#            if counts[b] > maxc:
#              maxb=b
#              maxc=counts[b]
#        graph(ys, Years, eval(maxb),fn1+"_"+str(iteration)+'.csv'+" "+str(maxc), picfn)
#        plt.close("all")
#
###==============================================================================
###NCDC Land, Ocean and Land/Ocean Hemispheric by SEASONS
###==============================================================================
  from NCDCset import NCDCFileSet
  fns=[
  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_LandOcean\\SH_LandOcean*1880-2015.csv',
  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_LandOcean\\NH_LandOcean*1880-2015.csv',
  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_LandOcean\\Global_LandOcean*1880-2015.csv',
  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Land\\SH-Land*1880-2015.csv',
  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Land\\NH-Land*1880-2015.csv',
  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Land\\Global-Land*1880-2015.csv',
  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Ocean\\SH-Ocean*1880-2015.csv',
  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Ocean\\NH-Ocean*1880-2015.csv',
  'C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\18mar2015_home\\NOAA_Ocean\\Global-Ocean*1880-2015.csv'
  ]
  for fn in fns:
    for SEASON in ["DJF","MAM","JJA","SON"]:    
      data = NCDCFileSet(fn)
      ys, YEARS =  data.seasonal("",SEASON)
      Years=np.array([int(d) for d in YEARS[:,1]],dtype=float)
      counts={}
      newdir= os.path.splitext((os.path.basename(data.filename(SEASON))))[0]
      mkdir_p(newdir)
      fn1=newdir+'//'+os.path.basename(data.filename(SEASON))
      picfn=os.path.basename(data.filename(SEASON))+"_pic"  
      if os.path.exists(picfn+".png"):
        print picfn+".png", "exists!!"
      else:
        for iteration in range(100):
          #data = NCDCFileSet(fn)
          xs = np.array([random.random() for y in ys])
          convergent_breaks.TraceFile = fn1+"_"+str(iteration)+'.trace'
          cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
          report(fn1+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [fn1])
          breaks=cb[1]
          sbreaks=str(breaks)
          if not sbreaks in counts:
            counts[sbreaks] = 0
          counts[sbreaks] += 1
        maxc=0
        maxb=''
        for b in counts.keys():
          if maxb=='':
            maxb=b
            maxc=counts[b]
          else:
            if counts[b] > maxc:
              maxb=b
              maxc=counts[b]
        graph(ys, Years, eval(maxb),fn1+"_"+str(iteration)+'.csv'+" "+str(maxc), picfn)
        plt.close("all")

#####################################################################################################
##The following block of code for NOAA's banded observation
#      
#  files = gatherNOAANames("C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\ZonalLandSeaASCII_NCDC")
#  #fn="C:\\Users\\s4493222\\Documents\\ReferenceData\\NOAA\\ZonalLandSeaASCII_NCDC\\aravg.mon.land_ocean.90S.90N.v3.5.4.201408.asc"
#  for fn in files:
#    newdir= os.path.splitext(os.path.basename(fn))[0]
#    fn1=newdir+"//"+os.path.basename(fn)
#    if os.path.exists(os.path.basename(fn1)+"_pic.png"):
#      print os.path.basename(fn1)+"_pic.png","exists!"
#    else:        
#      data=NOAAascii(fn).annually()
#      #print data.annually()
#      ys=np.array([row[2] for row in data]) #2014 not complete
#      Years=np.array([row[0] for row in data])
#      counts={}
#      mkdir_p(newdir)
#      for iteration in range(100):
#        fn2=fn1+"_"+str(iteration)
#        xs = np.array([random.random() for y in ys])
#        convergent_breaks.TraceFile = fn2+'.trace'
#        cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#        report(fn2+'.csv',cb, xs, ys, Years, [fn])
#        breaks=cb[1]
#        sbreaks=str(breaks)
#        if not sbreaks in counts:
#          counts[sbreaks] = 0
#        counts[sbreaks] += 1
#      maxc=0
#      maxb=''
#      for b in counts.keys():
#        if maxb=='':
#          maxb=b
#          maxc=counts[b]
#        else:
#          if counts[b] > maxc:
#            maxb=b
#            maxc=counts[b]
#      graph(ys, Years, eval(maxb),os.path.basename(fn1)+" "+str(maxc), os.path.basename(fn1)+"_pic")
#
##
########################################################################################################
####The following block of code for GISSTEMP3 banded observation
########################################################################################################
####The following block of code for GISSTEMP3 banded observation
#  c3=csvfile.CSVfile("C:\\Users\\s4493222\\Documents\\ReferenceData\\GISS\\GISSTEMPV3_Apr2015\\GISSTEMPv3-April2015.csv",skip_header=8)
#  for fn in c3.fields():
#    if fn != "Year":
#      fn1=os.path.basename(c3.filename())
#      if os.path.exists(fn1+"_"+fn+"_pic.png"):
#        print fn1+"_"+fn+"_pic.png","exists!"
#      else:
#        print fn1, fn
#        counts={}
#        newdir= os.path.splitext((fn1))[0]+"_"+fn
#        print newdir
#        mkdir_p(newdir)
#        ys=c3.data4name(fn)
#        Years=c3.time4name(fn)
#        for iteration in range(100):
#          xs = np.array([random.random() for y in ys])
#          convergent_breaks.TraceFile = newdir+"//"+fn1+"_"+fn+"_"+str(iteration)+'.trace'
#          cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#          report(newdir+"//"+fn1+"_"+fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [fn1,fn])
#          breaks=cb[1]
#          sbreaks=str(breaks)
#          if not sbreaks in counts:
#            counts[sbreaks] = 0
#          counts[sbreaks] += 1
#        maxc=0
#        maxb=''
#        for b in counts.keys():
#          if maxb=='':
#            maxb=b
#            maxc=counts[b]
#          else:
#            if counts[b] > maxc:
#              maxb=b
#              maxc=counts[b]
#        graph(ys, Years, eval(maxb),fn1+"_"+fn+"_"+str(99)+'.csv'+" "+str(maxc), fn1+"_"+fn+"_pic")
##older versions
#  c3=csvfile.CSVfile("GISSTEMPto6-2013b.csv",skip_header=9)
#  for fn in c3.fields():
#    if fn != "Year":
#      if os.path.exists(c3.filename()+"_"+fn+"_pic.png"):
#        print c3.filename()+"_"+fn+"_pic.png","exists!"
#      else:
#        counts={}
#        newdir= os.path.splitext((os.path.basename(c3.filename())))[0]+"_"+fn
#        print newdir
#        mkdir_p(newdir)
#        ys=c3.data4name(fn)
#        Years=c3.time4name(fn)
#        for iteration in range(100):
#          xs = np.array([random.random() for y in ys])
#          convergent_breaks.TraceFile = newdir+"//"+c3.filename()+"_"+fn+"_"+str(iteration)+'.trace'
#          cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#          report(newdir+"//"+c3.filename()+"_"+fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [c3.filename(),fn])
#          breaks=cb[1]
#          sbreaks=str(breaks)
#          if not sbreaks in counts:
#            counts[sbreaks] = 0
#          counts[sbreaks] += 1
#        maxc=0
#        maxb=''
#        for b in counts.keys():
#          if maxb=='':
#            maxb=b
#            maxc=counts[b]
#          else:
#            if counts[b] > maxc:
#              maxb=b
#              maxc=counts[b]
#        graph(ys, Years, eval(maxb),c3.filename()+"_"+fn+"_"+str(99)+'.csv'+" "+str(maxc), c3.filename()+"_"+fn+"_pic")
##and the difference between the two sets
#  fngiss="C:\\Users\\s4493222\\Documents\\ReferenceData\\GISS\\GISSTEMPV3_Apr2015\\GISSTEMP_v3_April_2015_withDiffs_from_2014.csv"
#  c3=csvfile.CSVfile(fngiss,skip_header=0)
#  for fn in c3.fields():
#    if fn != "Year":
#      fn1=os.path.basename(c3.filename())
#      if os.path.exists(fn1+"_"+fn+"_pic.png"):
#        print fn1+"_"+fn+"_pic.png","exists!"
#      else:
#        print fn1, fn
#        counts={}
#        newdir= os.path.splitext((fn1))[0]+"_"+fn
#        print newdir
#        mkdir_p(newdir)
#        ys=c3.data4name(fn)
#        Years=c3.time4name(fn)
#        for iteration in range(100):
#          xs = np.array([random.random() for y in ys])
#          convergent_breaks.TraceFile = newdir+"//"+fn1+"_"+fn+"_"+str(iteration)+'.trace'
#          cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#          report(newdir+"//"+fn1+"_"+fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [fn1,fn])
#          breaks=cb[1]
#          sbreaks=str(breaks)
#          if not sbreaks in counts:
#            counts[sbreaks] = 0
#          counts[sbreaks] += 1
#        maxc=0
#        maxb=''
#        for b in counts.keys():
#          if maxb=='':
#            maxb=b
#            maxc=counts[b]
#          else:
#            if counts[b] > maxc:
#              maxb=b
#              maxc=counts[b]
#        graph(ys, Years, eval(maxb),fn1+"_"+fn+"_"+str(99)+'.csv'+" "+str(maxc), fn1+"_"+fn+"_pic")
#
#        
##
###  collated(".\\",0,initial_header=9)      
################################################################################
##The following code for HADCRUT4
##For 4.2
#  files=glob.glob('C:\\\Users\\s4493222\\Documents\\ReferenceData\\Hadley_Data\HADCRUT4\\HadCRUT.4.2.0.0.annual*.txt')
#  for fn1 in files:
#    if os.path.exists(os.path.basename(fn1)+"_pic.png"):
#      print os.path.basename(fn1)+"_pic.png","exists!"
#    else:
#        
#      counts={}
#      newdir= os.path.splitext((os.path.basename(fn1)))[0]
#      mkdir_p(newdir)
#      had=HADCRUT4(fn1)
#      Years=had.years()
#      ys=had.annual()
#      fn = newdir+"\\"+os.path.basename(had.filename())
#      for iteration in range(100):
#        xs = np.array([random.random() for y in ys])
#        convergent_breaks.TraceFile = fn+"_"+str(iteration)+'.trace'
#        cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#        report(fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [had.filename(),fn])
#        breaks=cb[1]
#        sbreaks=str(breaks)
#        if not sbreaks in counts:
#          counts[sbreaks] = 0
#        counts[sbreaks] += 1
#      maxc=0
#      maxb=''
#      for b in counts.keys():
#        if maxb=='':
#          maxb=b
#          maxc=counts[b]
#        else:
#          if counts[b] > maxc:
#            maxb=b
#            maxc=counts[b]
#      graph(ys, Years, eval(maxb),os.path.basename(fn)+" "+str(maxc), os.path.basename(fn)+"_pic")
#      plt.close("all")
#
##For 4.3.0.0
#  files=glob.glob('C:\\\Users\\s4493222\\Documents\\ReferenceData\\Hadley_Data\HADCRUT4\\HadCRUT.4.3.0.0.annual*.txt')
#  
#  for fn1 in files:
#    if os.path.exists(os.path.basename(fn1)+"_pic.png"):
#      print os.path.basename(fn1)+"_pic.png","exists!"
#    else:
#        
#      counts={}
#      newdir= os.path.splitext((os.path.basename(fn1)))[0]
#      mkdir_p(newdir)
#      had=HADCRUT4(fn1)
#      Years=had.years()
#      ys=had.annual()
#      fn = newdir+"\\"+os.path.basename(had.filename())
#      for iteration in range(100):
#        xs = np.array([random.random() for y in ys])
#        convergent_breaks.TraceFile = fn+"_"+str(iteration)+'.trace'
#        cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#        report(fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [had.filename(),fn])
#        breaks=cb[1]
#        sbreaks=str(breaks)
#        if not sbreaks in counts:
#          counts[sbreaks] = 0
#        counts[sbreaks] += 1
#      maxc=0
#      maxb=''
#      for b in counts.keys():
#        if maxb=='':
#          maxb=b
#          maxc=counts[b]
#        else:
#          if counts[b] > maxc:
#            maxb=b
#            maxc=counts[b]
#      graph(ys, Years, eval(maxb),os.path.basename(fn)+" "+str(maxc), os.path.basename(fn)+"_pic")
#      plt.close("all")
#########################################################################################
##This block of code for Hadsst
##For 3.1.1 
#  files=glob.glob('C:\\Users\\s4493222\\Documents\\ReferenceData\\Hadley_Data\\18Mar2015_home\\HadSST.3.1.1.0_*_ts.txt')
#  
#  for fn1 in files:
#    if os.path.exists(os.path.basename(fn1)+"_pic.png"):
#      print os.path.basename(fn1)+"_pic.png","exists!"
#    else:
#        
#      counts={}
#      newdir= os.path.splitext((os.path.basename(fn1)))[0]
#      mkdir_p(newdir)
#      had=HADCRUT4(fn1,annual=False)
#      Years=had.years()
#      ys=had.annual()
#      fn = newdir+"\\"+os.path.basename(had.filename())
#      for iteration in range(100):
#        xs = np.array([random.random() for y in ys])
#        convergent_breaks.TraceFile = fn+"_"+str(iteration)+'.trace'
#        cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#        report(fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [had.filename(),fn])
#        breaks=cb[1]
#        sbreaks=str(breaks)
#        if not sbreaks in counts:
#          counts[sbreaks] = 0
#        counts[sbreaks] += 1
#      maxc=0
#      maxb=''
#      for b in counts.keys():
#        if maxb=='':
#          maxb=b
#          maxc=counts[b]
#        else:
#          if counts[b] > maxc:
#            maxb=b
#            maxc=counts[b]
#      graph(ys, Years, eval(maxb),os.path.basename(fn)+" "+str(maxc), os.path.basename(fn)+"_pic")
#      plt.close("all")
###and monthly of the same TOBE COMPLETED
#  files=glob.glob('C:\\Users\\s4493222\\Documents\\ReferenceData\\Hadley_Data\\18Mar2015_home\\HadSST.3.1.1.0_*_ts.txt')
#  Labels=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
#  for month in range(12):      
#    for fn1 in files:
#      if os.path.exists(os.path.basename(fn1)+"_"+Labels[month]+"_pic.png"):
#        print os.path.basename(fn1)+"_"+Labels[month]+"_pic.png","exists!"
#      else:
#        print os.path.basename(fn1)+"_"+Labels[month]
#        counts={}
#        newdir= os.path.splitext((os.path.basename(fn1)))[0]+"_"+Labels[month]
#        mkdir_p(newdir)
#        had=HADCRUT4(fn1,annual=False)
#        Years=had.years()
#        ys=had.monthly(month)
#        fn = newdir+"\\"+os.path.basename(had.filename(Labels[month]))
#        for iteration in range(100):
#          xs = np.array([random.random() for y in ys])
#          convergent_breaks.TraceFile = fn+"_"+str(iteration)+'.trace'
#          cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#          report(fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [had.filename(Labels[month]),fn])
#          breaks=cb[1]
#          sbreaks=str(breaks)
#          if not sbreaks in counts:
#            counts[sbreaks] = 0
#          counts[sbreaks] += 1
#        maxc=0
#        maxb=''
#        for b in counts.keys():
#          if maxb=='':
#            maxb=b
#            maxc=counts[b]
#          else:
#            if counts[b] > maxc:
#              maxb=b
#              maxc=counts[b]
#        graph(ys, Years, eval(maxb),os.path.basename(fn)+" "+str(maxc), os.path.basename(fn)+"_pic")
#        plt.close("all")
        
##==============================================================================
##This code for the SEASONAL CRU compilations of HadSST, CRUTem and HadCRUT4.3
##==============================================================================
  files=glob.glob('C:\\Users\\s4493222\\Documents\\ReferenceData\\Hadley_Data\\18Mar2015_home\\FromCRU\\*.dat.txt')
  Labels=["DJF","MAM","JJA","SON"]
  for season in Labels:  
    for fn1 in files:
      if os.path.exists(os.path.basename(fn1)+"_"+season+"_pic.png"):
        print os.path.basename(fn1)+"_"+season+"_pic.png","exists!"
      else:
        counts={}
        newdir= os.path.splitext((os.path.basename(fn1)))[0]+"_"+season
        mkdir_p(newdir)
        had=CRU(fn1)
        Years=had.years()
        ys, dates=had.seasonal("",season)
        Years=np.array([int(d) for d in dates[:,1]],dtype=float)
        
        l=min(len(Years), len(ys))
        Years=Years[:l]
        ys=ys[:l]
        fn = newdir+"\\"+os.path.basename(had.filename())+"_"+season
        for iteration in range(100):
          xs = np.array([random.random() for y in ys])
          convergent_breaks.TraceFile = fn+"_"+str(iteration)+'.trace'
          cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
          report(fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [had.filename(),fn])
          breaks=cb[1]
          sbreaks=str(breaks)
          if not sbreaks in counts:
            counts[sbreaks] = 0
          counts[sbreaks] += 1
        maxc=0
        maxb=''
        for b in counts.keys():
          if maxb=='':
            maxb=b
            maxc=counts[b]
          else:
            if counts[b] > maxc:
              maxb=b
              maxc=counts[b]
        graph(ys, Years, eval(maxb),os.path.basename(fn)+" "+str(maxc), os.path.basename(fn)+"_pic")
        plt.close("all")
##==============================================================================
##This code for the monthly and annual CRU compilations of HadSST, CRUTem and HadCRUT4.3
##==============================================================================
#  files=glob.glob('C:\\Users\\s4493222\\Documents\\ReferenceData\\Hadley_Data\\18Mar2015_home\\FromCRU\\*.dat.txt')
#  Labels=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
#  for fn1 in files:
#    if os.path.exists(os.path.basename(fn1)+"_pic.png"):
#      print os.path.basename(fn1)+"_pic.png","exists!"
#    else:
#        
#      counts={}
#      newdir= os.path.splitext((os.path.basename(fn1)))[0]
#      mkdir_p(newdir)
#      had=CRU(fn1)
#      Years=had.years()
#      ys=had.annual()
#      fn = newdir+"\\"+os.path.basename(had.filename())
#      for iteration in range(100):
#        xs = np.array([random.random() for y in ys])
#        convergent_breaks.TraceFile = fn+"_"+str(iteration)+'.trace'
#        cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#        report(fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [had.filename(),fn])
#        breaks=cb[1]
#        sbreaks=str(breaks)
#        if not sbreaks in counts:
#          counts[sbreaks] = 0
#        counts[sbreaks] += 1
#      maxc=0
#      maxb=''
#      for b in counts.keys():
#        if maxb=='':
#          maxb=b
#          maxc=counts[b]
#        else:
#          if counts[b] > maxc:
#            maxb=b
#            maxc=counts[b]
#      graph(ys, Years, eval(maxb),os.path.basename(fn)+" "+str(maxc), os.path.basename(fn)+"_pic")
#      plt.close("all")
#  for month in range(12):      
#    for fn1 in files:
#      if os.path.exists(os.path.basename(fn1)+"_"+Labels[month]+"_pic.png"):
#        print os.path.basename(fn1)+"_"+Labels[month]+"_pic.png","exists!"
#      else:
#        print os.path.basename(fn1)+"_"+Labels[month]
#        counts={}
#        newdir= os.path.splitext((os.path.basename(fn1)))[0]+"_"+Labels[month]
#        mkdir_p(newdir)
#        had=CRU(fn1)
#        Years=had.years()
#        ys=had.monthly(month)
#        fn = newdir+"\\"+os.path.basename(had.filename(Labels[month]))
#        for iteration in range(100):
#          xs = np.array([random.random() for y in ys])
#          convergent_breaks.TraceFile = fn+"_"+str(iteration)+'.trace'
#          cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#          report(fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [had.filename(Labels[month]),fn])
#          breaks=cb[1]
#          sbreaks=str(breaks)
#          if not sbreaks in counts:
#            counts[sbreaks] = 0
#          counts[sbreaks] += 1
#        maxc=0
#        maxb=''
#        for b in counts.keys():
#          if maxb=='':
#            maxb=b
#            maxc=counts[b]
#          else:
#            if counts[b] > maxc:
#              maxb=b
#              maxc=counts[b]
#        graph(ys, Years, eval(maxb),os.path.basename(fn)+" "+str(maxc), os.path.basename(fn)+"_pic")
#        plt.close("all")
#        
###==============================================================================
###Cowtan and Way
###==============================================================================
#  from COWTAN_WAY import COWTAN_WAY
#  fn='C:\Users\s4493222\Documents\ReferenceData\CowtanWay\coverage2013\had4_krig_annual_v2_0_0.txt'
#  data = COWTAN_WAY(fn)
#  ys =  data.annual()
#  Years= data.years()
#  counts={}
#  newdir= os.path.splitext((os.path.basename(data.filename())))[0]
#  mkdir_p(newdir)
#  fn1=newdir+'//'+os.path.basename(data.filename())
#  picfn=os.path.basename(data.filename())+"_pic"  
#  if os.path.exists(picfn+".png"):
#    print picfn+".png", "exists!!"
#  else:
#    for iteration in range(100):
#      data = COWTAN_WAY(fn)
#      ys =  data.annual()
#      xs = np.array([random.random() for y in ys])
#      convergent_breaks.TraceFile = fn1+"_"+str(iteration)+'.trace'
#      cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#      report(fn1+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [fn1])
#      breaks=cb[1]
#      sbreaks=str(breaks)
#      if not sbreaks in counts:
#        counts[sbreaks] = 0
#      counts[sbreaks] += 1
#    maxc=0
#    maxb=''
#    for b in counts.keys():
#      if maxb=='':
#        maxb=b
#        maxc=counts[b]
#      else:
#        if counts[b] > maxc:
#          maxb=b
#          maxc=counts[b]
#    graph(ys, Years, eval(maxb),fn1+"_"+str(iteration)+'.csv'+" "+str(maxc), picfn)
#    plt.close("all")
#    
###    
####==============================================================================
####BEST
####==============================================================================
###    
#  import BERKLEY
#  data = BERKLEY.BERKLEY()
#  fn=data.filename()
#  ys =  data.annual()
#  Years= data.years()
#  counts={}
#  newdir= os.path.splitext((os.path.basename(data.filename())))[0]
#  mkdir_p(newdir)
#  fn1=newdir+'//'+os.path.basename(data.filename())
#  picfn=os.path.basename(data.filename())+"_pic"  
#  if os.path.exists(picfn+".png"):
#    print picfn+".png", "exists!!"
#  else:
#    for iteration in range(100):
#      xs = np.array([random.random() for y in ys])
#      convergent_breaks.TraceFile = fn1+"_"+str(iteration)+'.trace'
#      cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#      report(fn1+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [fn1])
#      breaks=cb[1]
#      sbreaks=str(breaks)
#      if not sbreaks in counts:
#        counts[sbreaks] = 0
#      counts[sbreaks] += 1
#    maxc=0
#    maxb=''
#    for b in counts.keys():
#      if maxb=='':
#        maxb=b
#        maxc=counts[b]
#      else:
#        if counts[b] > maxc:
#          maxb=b
#          maxc=counts[b]
#    graph(ys, Years, eval(maxb),fn1+"_"+str(iteration)+'.csv'+" "+str(maxc), picfn)
#    plt.close("all")
###    
####==============================================================================
####FOSTER RAHMSTORF 2011
####==============================================================================
#  fnin=os.environ["HOMEPATH"]+"\\Documents\\ReferenceData\\foster_rahmstorf_2011\\adjusted.csv"  
#  c3=csvfile.CSVfile(fnin)
#  for field in c3.fields():
#    if field != "t":
#      if TRENDS:
#        fn=fnin+"_T2_SCREEN_"+field
#      else:
#        fn=fnin+"_"+field
#      ys =  c3[field]
#      Years= c3["t"]
#      Years=np.array([int(np.mean(Years[i * 12: (i+1)*12])) for i in range(int(len(ys)/12))])
#      ys=np.array([np.mean(ys[i * 12: (i+1)*12]) for i in range(len(Years))])
#      counts={}
#      newdir= os.path.basename(fn)
#      mkdir_p(newdir)
#      fn1=newdir+'//'+os.path.basename(fn)
#      picfn=os.path.basename(fn)+"_pic"  
#      if os.path.exists(picfn+".png"):
#        print picfn+".png", "exists!!"
#      else:
#        for iteration in range(100):
#          xs = np.array([random.random() for y in ys])
#          convergent_breaks.TraceFile = fn1+"_"+str(iteration)+'.trace'
#          #cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#          cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#          report(fn1+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [fnin])
#          breaks=cb[1]
#          sbreaks=str(breaks)
#          if not sbreaks in counts:
#            counts[sbreaks] = 0
#          counts[sbreaks] += 1
#        maxc=0
#        maxb=''
#        for b in counts.keys():
#          if maxb=='':
#            maxb=b
#            maxc=counts[b]
#          else:
#            if counts[b] > maxc:
#              maxb=b
#              maxc=counts[b]
#        graph(ys, Years, eval(maxb),os.path.basename(fn) +" "+str(maxc), picfn)
#        plt.close("all")
##  
###==============================================================================
###Roger's test set
###==============================================================================
##    
#  fnin=os.environ["HOMEPATH"]+"\\Documents\\abrupt\\CMIP5_breakpoints\\artifishul\\Artyfishul_dartar4Jim_(1).csv"
#  c3=csvfile.CSVfile(fnin)
#  for field in c3.fields():
#    if field != "Year":
#      if TRENDS:
#        fn=fnin+"_T2_SCREEN_"+field
#      else:
#        fn=fnin+"_"+field
#      ys =  c3[field]
#      Years= c3["Year"]
#      counts={}
#      newdir= os.path.basename(fn)
#      mkdir_p(newdir)
#      fn1=newdir+'//'+os.path.basename(fn)
#      picfn=os.path.basename(fn)+"_pic"  
#      if os.path.exists(picfn+".png"):
#        print picfn+".png", "exists!!"
#      else:
#        for iteration in range(100):
#          xs = np.array([random.random() for y in ys])
#          convergent_breaks.TraceFile = fn1+"_"+str(iteration)+'.trace'
##          cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#          cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#          report(fn1+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [fnin])
#          breaks=cb[1]
#          sbreaks=str(breaks)
#          if not sbreaks in counts:
#            counts[sbreaks] = 0
#          counts[sbreaks] += 1
#        maxc=0
#        maxb=''
#        for b in counts.keys():
#          if maxb=='':
#            maxb=b
#            maxc=counts[b]
#          else:
#            if counts[b] > maxc:
#              maxb=b
#              maxc=counts[b]
#        graph(ys, Years, eval(maxb),os.path.basename(fn) +" "+str(maxc), picfn)
#        plt.close("all")
#
#
#
######################################################################################################
#
##The following block of code for ICMP data set sent by Roger
##controls
#  files=glob.glob('C:\\\Users\\s4493222\\Documents\\ReferenceData\\icmip5\\r*_piCont*.dat')
#  for fn1 in files:
#    if os.path.exists(os.path.basename(fn1)+"_pic.png"):
#      print os.path.basename(fn1)+"_pic.png","exists!"
#    else:
#        
#      counts={}
#      newdir= os.path.splitext((os.path.basename(fn1)))[0]
#      mkdir_p(newdir)
#      icmp=ICMP(fn1)
#      Years=icmp.years()
#      ys=icmp.annual_means()
#      fn = newdir+"\\"+os.path.basename(icmp.filename())
#      for iteration in range(100):
#        xs = np.array([random.random() for y in ys])
#        convergent_breaks.TraceFile = fn+"_"+str(iteration)+'.trace'
#        cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#        report(fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [icmp.filename(),fn])
#        breaks=cb[1]
#        sbreaks=str(breaks)
#        if not sbreaks in counts:
#          counts[sbreaks] = 0
#        counts[sbreaks] += 1
#      maxc=0
#      maxb=''
#      for b in counts.keys():
#        if maxb=='':
#          maxb=b
#          maxc=counts[b]
#        else:
#          if counts[b] > maxc:
#            maxb=b
#            maxc=counts[b]
#      graph(ys, Years, eval(maxb),os.path.basename(fn)+" "+str(maxc), os.path.basename(fn)+"_pic")
#      plt.close("all")
#  #then all RCPs 
#  files=glob.glob('C:\\\Users\\s4493222\\Documents\\ReferenceData\\icmip5\\r*_RCP*.dat')
#  for fn1 in files:
#    if os.path.exists(os.path.basename(fn1)+"_pic.png"):
#      print os.path.basename(fn1)+"_pic.png","exists!"
#    else:
#        
#      counts={}
#      newdir= os.path.splitext((os.path.basename(fn1)))[0]
#      mkdir_p(newdir)
#      icmp=ICMP(fn1)
#      Years=icmp.years()
#      ys=icmp.annual_means()
#      fn = newdir+"\\"+os.path.basename(icmp.filename())
#      for iteration in range(100):
#        xs = np.array([random.random() for y in ys])
#        convergent_breaks.TraceFile = fn+"_"+str(iteration)+'.trace'
#        cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#        report(fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [icmp.filename(),fn])
#        breaks=cb[1]
#        sbreaks=str(breaks)
#        if not sbreaks in counts:
#          counts[sbreaks] = 0
#        counts[sbreaks] += 1
#      maxc=0
#      maxb=''
#      for b in counts.keys():
#        if maxb=='':
#          maxb=b
#          maxc=counts[b]
#        else:
#          if counts[b] > maxc:
#            maxb=b
#            maxc=counts[b]
#      graph(ys, Years, eval(maxb),os.path.basename(fn)+" "+str(maxc), os.path.basename(fn)+"_pic")
#      plt.close("all")
#      
#####################################################################################################
##The following block of code for KNMI data downloaded by jim to match set sent by Roger
#
#  files=glob.glob('C:\\\Users\\s4493222\\Documents\\ReferenceData\\icmip5\\rcp26\\r*.dat.txt')
#  for fn1 in files:
#    if os.path.exists(os.path.basename(fn1)+"_pic.png"):
#      print os.path.basename(fn1)+"_pic.png","exists!"
#    else:
#        
#      counts={}
#      newdir= os.path.splitext((os.path.basename(fn1)))[0]
#      mkdir_p(newdir)
#      icmp=ICMP(fn1)
#      Years=icmp.years()
#      ys=icmp.annual_means()
#      fn = newdir+"\\"+os.path.basename(icmp.filename())
#      for iteration in range(100):
#        xs = np.array([random.random() for y in ys])
#        convergent_breaks.TraceFile = fn+"_"+str(iteration)+'.trace'
#        cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#        report(fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [icmp.filename(),fn])
#        breaks=cb[1]
#        sbreaks=str(breaks)
#        if not sbreaks in counts:
#          counts[sbreaks] = 0
#        counts[sbreaks] += 1
#      maxc=0
#      maxb=''
#      for b in counts.keys():
#        if maxb=='':
#          maxb=b
#          maxc=counts[b]
#        else:
#          if counts[b] > maxc:
#            maxb=b
#            maxc=counts[b]
#      graph(ys, Years, eval(maxb),os.path.basename(fn)+" "+str(maxc), os.path.basename(fn)+"_pic")
#      plt.close("all")
#      
###==============================================================================
###cmip 3 files from QCCCE
###==============================================================================
##    
#  files =gatherCMIP3Names("C:\\Users\\s4493222\\Documents\\ReferenceData\\CSIRO_Model_Set\\gwFiles")
#  for fnx in files:
#    fn1=fnx
#    if os.path.exists(os.path.basename(fn1)+"_pic.png"):
#      print os.path.basename(fn1)+"_pic.png","exists!"
#    else:
#      print fnx
#      counts={}
#      newdir= os.path.splitext((os.path.basename(fn1)))[0]
#      mkdir_p(newdir)
#      c3=CMIP3gw(fnx)
#      ys =  c3.Warming()
#      xs = np.array([random.random() for y in ys])
#      Years= c3.Years()
#      fn = newdir+"\\"+os.path.basename(fn1)
#      for iteration in range(100):
#        xs = np.array([random.random() for y in ys])
#        convergent_breaks.TraceFile = fn+"_"+str(iteration)+'.trace'
#        cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#        report(fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [fn1,fn])
#        breaks=cb[1]
#        sbreaks=str(breaks)
#        if not sbreaks in counts:
#          counts[sbreaks] = 0
#        counts[sbreaks] += 1
#      maxc=0
#      maxb=''
#      for b in counts.keys():
#        if maxb=='':
#          maxb=b
#          maxc=counts[b]
#        else:
#          if counts[b] > maxc:
#            maxb=b
#            maxc=counts[b]
#      graph(ys, Years, eval(maxb),os.path.basename(fn)+" "+str(maxc), os.path.basename(fn)+"_pic")
#      plt.close("all")
#  '''
#  Traceback (most recent call last):
#    File "<stdin>", line 1, in <module>
#    File "C:\Users\Public\Documents\bin\WinPython-32bit-2.7.6.2\python-2.7.6\lib\site-packages\spyderlib\widgets\externalshell\sitecustomize.py", line 540, in runfile
#      execfile(filename, namespace)
#    File "C:/Users/s4493222/Documents/abrupt/4Roger_ClimatePapers_SVN_280/CB_Reports.py", line 1360, in <module>
#      cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#    File "ConvergentBreaks.py", line 614, in convergentBreaks
#      return convergentBreaks_Inner(testdata, controldata, datayears, controldata, model, pr=pr, screenpr=screenpr, trace=trace, shallow=shallow)
#    File "ConvergentBreaks.py", line 542, in convergentBreaks_Inner
#      tf=open(TraceFile,"w")
#  IOError: [Errno 2] No such file or directory: 'ncar_ccsm3_0.sresa1b.run2.monthly.tas_A1.SRESA1B_2.CCSM.atmm.2000-01_cat_2099-12_1870-2099_gw\\ncar_ccsm3_0.sresa1b.run2.monthly.tas_A1.SRESA1B_2.CCSM.atmm.2000-01_cat_2099-12_1870-2099_gw.txt_0.trace'
#  '''
###==============================================================================
###cmip 5 files from QCCCE
###==============================================================================
##    
##
#  files=gatherCMIP5Names('C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW')
#  for fnx in files:
#    fn1=fnx[1]
#    if os.path.exists(os.path.basename(fn1)+"_pic.png"):
#      print os.path.basename(fn1)+"_pic.png","exists!"
#    else:
#      print fnx
#      counts={}
#      newdir= os.path.splitext((os.path.basename(fn1)))[0]
#      mkdir_p(newdir)
#      c3=CMIP5gw(fnx)
#      ys =  c3.Warming()
#      xs = np.array([random.random() for y in ys])
#      Years= c3.Years()
#      fn = newdir+"\\"+os.path.basename(fn1)
#      for iteration in range(100):
#        xs = np.array([random.random() for y in ys])
#        convergent_breaks.TraceFile = fn+"_"+str(iteration)+'.trace'
#        cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#        report(fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [fn1,fn])
#        breaks=cb[1]
#        sbreaks=str(breaks)
#        if not sbreaks in counts:
#          counts[sbreaks] = 0
#        counts[sbreaks] += 1
#      maxc=0
#      maxb=''
#      for b in counts.keys():
#        if maxb=='':
#          maxb=b
#          maxc=counts[b]
#        else:
#          if counts[b] > maxc:
#            maxb=b
#            maxc=counts[b]
#      graph(ys, Years, eval(maxb),os.path.basename(fn)+" "+str(maxc), os.path.basename(fn)+"_pic")
#      plt.close("all")
   
#######################################################################################################

#==============================================================================
#   files =gatherCMIP3Names("C:\\Users\\s4493222\\Documents\\ReferenceData\\CSIRO_Model_Set\\gwFiles")
#   badf="C:\\Users\\s4493222\\Documents\\ReferenceData\\CSIRO_Model_Set\\gwFiles\\cnrm_cm3.sresa2.run1.monthly.tas_A1_1860-2099_gw.txt"
#   for fn in files:#[badf]:#files:
#     c3=CMIP3gw(fn)
#     ys =  c3.Warming()
#     xs = np.array([random.random() for y in ys])
#     Years= c3.Years()
#     try:
#       convergent_breaks.TraceFile = fn+'.trace'
#       cb=convergent_breaks.convergentBreaks(ys, xs, Years, os.path.basename(fn), mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#       
#     except Exception as e:
#       #convergent_breaks.tf.close()
#       print  os.path.basename(fn)
#       raise
#     report(fn+'.csv',cb, xs, ys, Years, [fn])


#  
#  DATA=np.array([[0 ,-0.399 ,0.604823921142 ,1850.0 ,0.604823921142],
#  [1 ,-0.219 ,0.767573180529 ,1851.0 ,0.767573180529],
#  [2 ,-0.243 ,0.190444637065 ,1852.0 ,0.190444637065],
#  [3 ,-0.313 ,0.300456131497 ,1853.0 ,0.300456131497],
#  [4 ,-0.235 ,0.0491362919722 ,1854.0 ,0.0491362919722],
#  [5 ,-0.233 ,0.763290296809 ,1855.0 ,0.763290296809],
#  [6 ,-0.39 ,0.882352600577 ,1856.0 ,0.882352600577],
#  [7 ,-0.467 ,0.346831189493 ,1857.0 ,0.346831189493],
#  [8 ,-0.42 ,0.41158011085 ,1858.0 ,0.41158011085],
#  [9 ,-0.369 ,0.845582916861 ,1859.0 ,0.845582916861],
#  [10 ,-0.452 ,0.850988506537 ,1860.0 ,0.850988506537],
#  [11 ,-0.488 ,0.656150822511 ,1861.0 ,0.656150822511],
#  [12 ,-0.583 ,0.347359892438 ,1862.0 ,0.347359892438],
#  [13 ,-0.393 ,0.313595364009 ,1863.0 ,0.313595364009],
#  [14 ,-0.462 ,0.54041958302 ,1864.0 ,0.54041958302],
#  [15 ,-0.337 ,0.262452603112 ,1865.0 ,0.262452603112],
#  [16 ,-0.295 ,0.805539810623 ,1866.0 ,0.805539810623],
#  [17 ,-0.349 ,0.522741058078 ,1867.0 ,0.522741058078],
#  [18 ,-0.246 ,0.70477182259 ,1868.0 ,0.70477182259],
#  [19 ,-0.274 ,0.45796944645 ,1869.0 ,0.45796944645],
#  [20 ,-0.331 ,0.0156695949516 ,1870.0 ,0.0156695949516],
#  [21 ,-0.349 ,0.640499420097 ,1871.0 ,0.640499420097],
#  [22 ,-0.312 ,0.641413237403 ,1872.0 ,0.641413237403],
#  [23 ,-0.289 ,0.34896409242 ,1873.0 ,0.34896409242],
#  [24 ,-0.407 ,0.107375792302 ,1874.0 ,0.107375792302],
#  [25 ,-0.408 ,0.991185863843 ,1875.0 ,0.991185863843],
#  [26 ,-0.375 ,0.539191502725 ,1876.0 ,0.539191502725],
#  [27 ,-0.056 ,0.205798581487 ,1877.0 ,0.205798581487],
#  [28 ,0.029 ,0.308332905193 ,1878.0 ,0.308332905193],
#  [29 ,-0.263 ,0.872117971873 ,1879.0 ,0.872117971873],
#  [30 ,-0.274 ,0.197746939267 ,1880.0 ,0.197746939267],
#  [31 ,-0.229 ,0.972060123327 ,1881.0 ,0.972060123327],
#  [32 ,-0.277 ,0.342021697779 ,1882.0 ,0.342021697779],
#  [33 ,-0.322 ,0.798822597428 ,1883.0 ,0.798822597428],
#  [34 ,-0.483 ,0.211114290582 ,1884.0 ,0.211114290582],
#  [35 ,-0.454 ,0.126802905319 ,1885.0 ,0.126802905319],
#  [36 ,-0.443 ,0.449591861419 ,1886.0 ,0.449591861419],
#  [37 ,-0.482 ,0.970610755956 ,1887.0 ,0.970610755956],
#  [38 ,-0.346 ,0.096890237844 ,1888.0 ,0.096890237844],
#  [39 ,-0.207 ,0.163168167083 ,1889.0 ,0.163168167083],
#  [40 ,-0.454 ,0.44583181847 ,1890.0 ,0.44583181847],
#  [41 ,-0.381 ,0.481130550735 ,1891.0 ,0.481130550735],
#  [42 ,-0.469 ,0.749643183849 ,1892.0 ,0.749643183849],
#  [43 ,-0.444 ,0.078858184099 ,1893.0 ,0.078858184099],
#  [44 ,-0.413 ,0.0398562737092 ,1894.0 ,0.0398562737092],
#  [45 ,-0.381 ,0.580273258731 ,1895.0 ,0.580273258731],
#  [46 ,-0.217 ,0.479470609839 ,1896.0 ,0.479470609839],
#  [47 ,-0.226 ,0.453247404264 ,1897.0 ,0.453247404264],
#  [48 ,-0.446 ,0.0182359203687 ,1898.0 ,0.0182359203687],
#  [49 ,-0.292 ,0.664696172306 ,1899.0 ,0.664696172306],
#  [50 ,-0.183 ,0.239125612544 ,1900.0 ,0.239125612544],
#  [51 ,-0.237 ,0.0386045050558 ,1901.0 ,0.0386045050558],
#  [52 ,-0.397 ,0.438439950386 ,1902.0 ,0.438439950386],
#  [53 ,-0.468 ,0.698752292157 ,1903.0 ,0.698752292157],
#  [54 ,-0.524 ,0.358096389692 ,1904.0 ,0.358096389692],
#  [55 ,-0.36 ,0.949133149641 ,1905.0 ,0.949133149641],
#  [56 ,-0.273 ,0.205042495238 ,1906.0 ,0.205042495238],
#  [57 ,-0.442 ,0.334915046509 ,1907.0 ,0.334915046509],
#  [58 ,-0.471 ,0.90857431311 ,1908.0 ,0.90857431311],
#  [59 ,-0.522 ,0.711115494867 ,1909.0 ,0.711115494867],
#  [60 ,-0.505 ,0.574388040353 ,1910.0 ,0.574388040353],
#  [61 ,-0.526 ,0.989123422847 ,1911.0 ,0.989123422847],
#  [62 ,-0.454 ,0.862268927689 ,1912.0 ,0.862268927689],
#  [63 ,-0.423 ,0.819432226112 ,1913.0 ,0.819432226112],
#  [64 ,-0.247 ,0.366578721457 ,1914.0 ,0.366578721457],
#  [65 ,-0.171 ,0.31374345957 ,1915.0 ,0.31374345957],
#  [66 ,-0.404 ,0.875830628787 ,1916.0 ,0.875830628787],
#  [67 ,-0.495 ,0.0359072802989 ,1917.0 ,0.0359072802989],
#  [68 ,-0.368 ,0.00229800960199 ,1918.0 ,0.00229800960199],
#  [69 ,-0.308 ,0.801444102726 ,1919.0 ,0.801444102726],
#  [70 ,-0.268 ,0.419734309857 ,1920.0 ,0.419734309857],
#  [71 ,-0.207 ,0.905310679889 ,1921.0 ,0.905310679889],
#  [72 ,-0.293 ,0.206109868612 ,1922.0 ,0.206109868612],
#  [73 ,-0.269 ,0.722005415425 ,1923.0 ,0.722005415425],
#  [74 ,-0.263 ,0.83664367499 ,1924.0 ,0.83664367499],
#  [75 ,-0.23 ,0.783512324749 ,1925.0 ,0.783512324749],
#  [76 ,-0.088 ,0.128525634919 ,1926.0 ,0.128525634919],
#  [77 ,-0.194 ,0.368726318291 ,1927.0 ,0.368726318291],
#  [78 ,-0.167 ,0.159306071358 ,1928.0 ,0.159306071358],
#  [79 ,-0.346 ,0.32600988163 ,1929.0 ,0.32600988163],
#  [80 ,-0.131 ,0.0572803480239 ,1930.0 ,0.0572803480239],
#  [81 ,-0.071 ,0.34878232493 ,1931.0 ,0.34878232493],
#  [82 ,-0.113 ,0.0475059641984 ,1932.0 ,0.0475059641984],
#  [83 ,-0.279 ,0.809384115826 ,1933.0 ,0.809384115826],
#  [84 ,-0.135 ,0.457227724657 ,1934.0 ,0.457227724657],
#  [85 ,-0.172 ,0.0760701850103 ,1935.0 ,0.0760701850103],
#  [86 ,-0.128 ,0.14103805707 ,1936.0 ,0.14103805707],
#  [87 ,0.015 ,0.573902153851 ,1937.0 ,0.573902153851],
#  [88 ,0.012 ,0.846026625155 ,1938.0 ,0.846026625155],
#  [89 ,-0.013 ,0.571481003961 ,1939.0 ,0.571481003961],
#  [90 ,0.072 ,0.816056666008 ,1940.0 ,0.816056666008],
#  [91 ,0.059 ,0.0817146173678 ,1941.0 ,0.0817146173678],
#  [92 ,0.008 ,0.134858101998 ,1942.0 ,0.134858101998],
#  [93 ,0.045 ,0.615075265108 ,1943.0 ,0.615075265108],
#  [94 ,0.148 ,0.544670631024 ,1944.0 ,0.544670631024],
#  [95 ,0.02 ,0.569313775446 ,1945.0 ,0.569313775446],
#  [96 ,-0.066 ,0.171392723028 ,1946.0 ,0.171392723028],
#  [97 ,0.002 ,0.226766993276 ,1947.0 ,0.226766993276],
#  [98 ,-0.071 ,0.372147769762 ,1948.0 ,0.372147769762],
#  [99 ,-0.122 ,0.776028669334 ,1949.0 ,0.776028669334],
#  [100 ,-0.207 ,0.818348776426 ,1950.0 ,0.818348776426],
#  [101 ,-0.028 ,0.237852305661 ,1951.0 ,0.237852305661],
#  [102 ,0.046 ,0.586193286016 ,1952.0 ,0.586193286016],
#  [103 ,0.109 ,0.0189256874725 ,1953.0 ,0.0189256874725],
#  [104 ,-0.093 ,0.63541769868 ,1954.0 ,0.63541769868],
#  [105 ,-0.159 ,0.596121798528 ,1955.0 ,0.596121798528],
#  [106 ,-0.218 ,0.300992124876 ,1956.0 ,0.300992124876],
#  [107 ,0.011 ,0.0498158284434 ,1957.0 ,0.0498158284434],
#  [108 ,0.032 ,0.0169136092644 ,1958.0 ,0.0169136092644],
#  [109 ,-0.0 ,0.847119180678 ,1959.0 ,0.847119180678],
#  [110 ,-0.05 ,0.0841195308773 ,1960.0 ,0.0841195308773],
#  [111 ,0.025 ,0.477726910698 ,1961.0 ,0.477726910698],
#  [112 ,0.002 ,0.980808404719 ,1962.0 ,0.980808404719],
#  [113 ,0.035 ,0.899763453878 ,1963.0 ,0.899763453878],
#  [114 ,-0.242 ,0.557487506503 ,1964.0 ,0.557487506503],
#  [115 ,-0.141 ,0.82395696852 ,1965.0 ,0.82395696852],
#  [116 ,-0.075 ,0.297154128765 ,1966.0 ,0.297154128765],
#  [117 ,-0.04 ,0.965193491462 ,1967.0 ,0.965193491462],
#  [118 ,-0.097 ,0.0486540517973 ,1968.0 ,0.0486540517973],
#  [119 ,0.044 ,0.35118746892 ,1969.0 ,0.35118746892],
#  [120 ,-0.016 ,0.435686266083 ,1970.0 ,0.435686266083],
#  [121 ,-0.153 ,0.04614293683 ,1971.0 ,0.04614293683],
#  [122 ,-0.073 ,0.24881418928 ,1972.0 ,0.24881418928],
#  [123 ,0.06 ,0.193243723299 ,1973.0 ,0.193243723299],
#  [124 ,-0.176 ,0.683150435697 ,1974.0 ,0.683150435697],
#  [125 ,-0.121 ,0.10858828495 ,1975.0 ,0.10858828495],
#  [126 ,-0.231 ,0.494531889186 ,1976.0 ,0.494531889186],
#  [127 ,0.059 ,0.256634206604 ,1977.0 ,0.256634206604],
#  [128 ,-0.053 ,0.916298081392 ,1978.0 ,0.916298081392],
#  [129 ,0.032 ,0.964292986926 ,1979.0 ,0.964292986926],
#  [130 ,0.143 ,0.201982113697 ,1980.0 ,0.201982113697],
#  [131 ,0.19 ,0.648664919555 ,1981.0 ,0.648664919555],
#  [132 ,-0.006 ,0.14590252913 ,1982.0 ,0.14590252913],
#  [133 ,0.177 ,0.178348312781 ,1983.0 ,0.178348312781],
#  [134 ,0.007 ,0.469148310777 ,1984.0 ,0.469148310777],
#  [135 ,-0.006 ,0.183654648656 ,1985.0 ,0.183654648656],
#  [136 ,0.038 ,0.862872587382 ,1986.0 ,0.862872587382],
#  [137 ,0.167 ,0.8933103339 ,1987.0 ,0.8933103339],
#  [138 ,0.224 ,0.657414788633 ,1988.0 ,0.657414788633],
#  [139 ,0.112 ,0.684621549087 ,1989.0 ,0.684621549087],
#  [140 ,0.297 ,0.14391955398 ,1990.0 ,0.14391955398],
#  [141 ,0.289 ,0.517843606575 ,1991.0 ,0.517843606575],
#  [142 ,0.093 ,0.897457655364 ,1992.0 ,0.897457655364],
#  [143 ,0.135 ,0.741461443981 ,1993.0 ,0.741461443981],
#  [144 ,0.205 ,0.269149016924 ,1994.0 ,0.269149016924],
#  [145 ,0.332 ,0.652133117 ,1995.0 ,0.652133117],
#  [146 ,0.228 ,0.842834034841 ,1996.0 ,0.842834034841],
#  [147 ,0.381 ,0.832371028049 ,1997.0 ,0.832371028049],
#  [148 ,0.533 ,0.68753599042 ,1998.0 ,0.68753599042],
#  [149 ,0.29 ,0.703235005545 ,1999.0 ,0.703235005545],
#  [150 ,0.302 ,0.932904849719 ,2000.0 ,0.932904849719],
#  [151 ,0.456 ,0.89107558689 ,2001.0 ,0.89107558689],
#  [152 ,0.524 ,0.560927377569 ,2002.0 ,0.560927377569],
#  [153 ,0.525 ,0.583019028719 ,2003.0 ,0.583019028719],
#  [154 ,0.443 ,0.176960808231 ,2004.0 ,0.176960808231],
#  [155 ,0.59 ,0.909966214618 ,2005.0 ,0.909966214618],
#  [156 ,0.541 ,0.772952647972 ,2006.0 ,0.772952647972],
#  [157 ,0.568 ,0.246604572849 ,2007.0 ,0.246604572849],
#  [158 ,0.431 ,0.299854535057 ,2008.0 ,0.299854535057],
#  [159 ,0.561 ,0.532654263038 ,2009.0 ,0.532654263038],
#  [160 ,0.634 ,0.914854283023 ,2010.0 ,0.914854283023],
#  [161 ,0.496 ,0.477331071355 ,2011.0 ,0.477331071355],
#  [162 ,0.518 ,0.195098445278 ,2012.0 ,0.195098445278],
#  [163 ,0.534 ,0.934274338782 ,2013.0 ,0.934274338782]])
#==============================================================================
#  import statbreaks
# 
#  import bivariate_multi as bivariate
#  import random
# #  import bivariate
#  DATA=statbreaks.brkrpt("C:\\Users\\s4493222\\Documents\\abrupt\\4Roger_Nature_SVN_264\\GISSTEMPto6-2013b_SHem\\GISSTEMPto6-2013b.csv_SHem_99.trace")
#  results=[]
#  starts=[]
#  ys=DATA.data()[:,1][89:]
#  xs=DATA.data()[:,2][89:]
#  Years=DATA.data()[:,3][89:]
#  print Years[[0,-1]]
#  for i in range(100):
#     bv=bivariate.bivariate(ys,xs, anomalise=False, pr=0.01)
#     
#     #cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#     results.append(bv.maxIndexTi())
#     starts.append(bv.allPoints(pr=0.01)[1])
#     print results[-1],starts[-1]
#     xs = np.array([random.random() for i in ys])
#  print Years[list(set(results))], results, len(starts)
#==============================================================================
  
#  convergent_breaks.TraceFile = 'case20.trace'
#  for i in range(100):
#    cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#    results.append(cb[1])
#    starts.append(cb[0])
#    xs = np.array([random.random() for i in ys])
#  print results
  
#r1i1p1_inmcm4_icmip5_tas_Amon_ens_piControl_0-360E_-90-90N_n_su_030.final.csv
#r1i1p2_GISS-E2-H_RCP4.5_icmip5_tas_Amon_ens_rcp45to85_0-360E_-90-90N_n_su_052.final.csv


##==============================================================================
##cmip 5 files from QCCCE used as a quick test
##==============================================================================
#    
#
#  files=gatherCMIP5Names('C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW')
#  flist=[
#  "tas_Amon_ACCESS1-0_rcp85_r1i1p1_200601-210012.GW",
#  "tas_Amon_ACCESS1-3_rcp85_r1i1p1_200601-210012.GW",
#  "tas_Amon_bcc-csm1-1_rcp60_r1i1p1_200601-209912.GW",
#  "tas_Amon_bcc-csm1-1_rcp85_r1i1p1_200601-209912.GW",
#  "tas_Amon_CSIRO-Mk3-6-0_rcp26_r7i1p1_200601-210012.GW",
#  "tas_Amon_CSIRO-Mk3-6-0_rcp45_r1i1p1_200601-210012.GW",
#  "tas_Amon_CSIRO-Mk3-6-0_rcp45_r2i1p1_200601-210012.GW"]  
#  for fnx in files:
#    fn1=fnx[1]
#    doMe=False
#    for fl in flist:
#      if fn1.find(fl) >-1:
#        doMe = True
#        print "Doing ",fn1, fl 
#    if not doMe:
#      print "Not to do ",fn1
#    else:
#  #for fn in [['C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_bcc-csm1-1-m_historical_r1i1p1_185001-201212.GW', 'C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_bcc-csm1-1-m_rcp26_r1i1p1_200601-210012.GW']]:
#  #for fn in [['C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_CSIRO-Mk3-6-0_historical_r4i1p1_185001-200512.GW', 'C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW\\tas_Amon_CSIRO-Mk3-6-0_rcp26_r4i1p1_200601-210012.GW']]:
##    c3=CMIP5gw(fn1)
##    ys =  c3.Warming()
##    xs = np.array([random.random() for y in ys])
##    Years= c3.Years()
##    convergent_breaks.TraceFile = fn[1]+'.trace'
##    cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
##    report(fn1[1]+'.csv',cb, xs, ys, Years, fn)
##  for intype in ["rcp26","rcp45","rcp60","rcp85"]:
##    collatedCMIP5("C:\\Users\\s4493222\\Documents\\ReferenceData\\QCCCE_CMIP5_GW",intype)      
#      print fnx
#      counts={}
#      newdir= os.path.splitext((os.path.basename(fn1)))[0]
#      mkdir_p(newdir)
#      c3=CMIP5gw(fnx)
#      ys =  c3.Warming()
#      xs = np.array([random.random() for y in ys])
#      Years= c3.Years()
#      fn = newdir+"\\"+os.path.basename(fn1)
#      for iteration in range(100):
#        xs = np.array([random.random() for y in ys])
#        convergent_breaks.TraceFile = fn+"_"+str(iteration)+'.trace'
#        cb=convergent_breaks.convergentBreaks(ys, xs, Years, 'GWAnom', mode="control", guide="Stability",screenpr=SCREENPR, pr=0.01, trace=True)
#        report(fn+"_"+str(iteration)+'.csv',cb, xs, ys, Years, [fn1,fn])
#        breaks=cb[1]
#        sbreaks=str(breaks)
#        if not sbreaks in counts:
#          counts[sbreaks] = 0
#        counts[sbreaks] += 1
#      maxc=0
#      maxb=''
#      for b in counts.keys():
#        if maxb=='':
#          maxb=b
#          maxc=counts[b]
#        else:
#          if counts[b] > maxc:
#            maxb=b
#            maxc=counts[b]
#      graph(ys, Years, eval(maxb),os.path.basename(fn)+" "+str(maxc), os.path.basename(fn)+"_pic")
#      plt.close("all")
#   

#tas_Amon_ACCESS1-0_rcp85_r1i1p1_200601-210012.GW_pic 97
#tas_Amon_ACCESS1-3_rcp85_r1i1p1_200601-210012.GW_pic 33
#tas_Amon_bcc-csm1-1_rcp60_r1i1p1_200601-209912.GW_pic 26
#tas_Amon_bcc-csm1-1_rcp85_r1i1p1_200601-209912.GW_pic 13
#tas_Amon_CSIRO-Mk3-6-0_rcp26_r7i1p1_200601-210012.GW_pic 100
#tas_Amon_CSIRO-Mk3-6-0_rcp45_r1i1p1_200601-210012.GW_pic 98
#tas_Amon_CSIRO-Mk3-6-0_rcp45_r2i1p1_200601-210012.GW_pic 22

flist=[
"tas_Amon_ACCESS1-0_rcp85_r1i1p1_200601-210012.GW",
"tas_Amon_ACCESS1-3_rcp85_r1i1p1_200601-210012.GW",
"tas_Amon_bcc-csm1-1_rcp60_r1i1p1_200601-209912.GW",
"tas_Amon_bcc-csm1-1_rcp85_r1i1p1_200601-209912.GW",
"tas_Amon_CSIRO-Mk3-6-0_rcp26_r7i1p1_200601-210012.GW",
"tas_Amon_CSIRO-Mk3-6-0_rcp45_r1i1p1_200601-210012.GW",
"tas_Amon_CSIRO-Mk3-6-0_rcp45_r2i1p1_200601-210012.GW"]