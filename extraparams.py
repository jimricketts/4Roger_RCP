# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 19:00:04 2015

@author: James
"""
import statbreaks
import doancova
import os
import glob
import cmip5_info
import cmip3_info

SVNRevision="$Revision: 374 $"



def formatprobs(p1,p2=None):
  if p1 <=0.001:
    f1="****"
  elif p1 <= 0.01:
    f1="*** "
  elif p1 <= 0.05:
    f1 = "**  "
  elif p1<= 0.1:
    f1 = "*   "
  else: f1 = "N/S "

  if p2 != None:
    if p2 <=0.001:
      f2="****"
    elif p2 <= 0.01:
      f2="*** "
    elif p2 <= 0.05:
      f2 = "**  "
    elif p2<= 0.1:
      f2 = "*   "
    else: f2 = "N/S "
  
    return "%06f,%s, %06f,%s" % (p1,f1,p2,f2)
  else:
    return "%06f,%s" % (p1,f1)

def oldstyle():
  outf=open("report_2.txt", "w")
  with open(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\x.x_08Jun",'r') as f:
    lines=f.readlines()
  for line in lines:
    #print line
    pos1=line.find('"')
    if pos1>=0:
      pos2=line.find('"',pos1+1)
      #print pos1,pos2
      fn=line[:pos1-1]
      #print fn+":"
      yset=eval(line[pos1:pos2+1])
      #print "YSET", yset
      _,yr,pre,post=line[pos2+1:].split()[:4]
     # print yr,pre,post
      dfn=fn[:-len('final.csv')-1]
      tfn=os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\"+dfn
      #print fn, tfn 
      if dfn[:len("Global")]=="Global" or dfn[0] == 'r':
        print "do",tfn
        try:
          direc=glob.glob(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\"+fn[:-len('final.csv')-1]+'/*.trace')[0]
          bp=statbreaks.brkrpt(direc)
          #print bp.years()
          yrix=int(list(bp.years()).index(float(yr)))
          preix=int(list(bp.years()).index(float(pre)))
          postix=int(list(bp.years()).index(float(post)))
          #print bp.years()[preix:postix+1]
          ancova=doancova.ANCOVA(bp.years()[preix:postix+1], bp.ys()[preix:postix+1], bp.years()[yrix+1],Rpath='"C:\\Program Files\\R\\R-3.0.2\\bin\\Rscript.exe"')
          ancova.parseResult(ancova.results())
          print >>outf, "slope",line[:-1],":",ancova._ANCOVA__breakYear,ancova.firstSlope(),formatprobs(ancova.firstSlopeProb()), ancova.firstRsq(), ancova.secondSlope(),formatprobs(ancova.secondSlopeProb()),ancova.secondRsq(),ancova.shift(),formatprobs(ancova._ANCOVA__prDiffShift),":",
          if dfn[:len("Global")]=="Global":
            print >>outf, cmip3_info.matchfn(tfn),
          elif dfn[0] == 'r':
            print >>outf, cmip5_info.matchfn(tfn),
          print >>outf, formatprobs(ancova.changeOfTrendProb(),ancova._ANCOVA__prDiffShift),formatprobs(ancova.changeOfTrendProb()*ancova._ANCOVA__prDiffShift)
          
          #if 
        except Exception as e:
          fl=os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\"+fn[:-len('final.csv')-1]+'/*.trace'
          print >>outf,"\n",line,"Error getting files",fl
          print "using ", fl
          #print glob.glob(fl)
          #raise Exception(str(e))
          pass
        
      #print direc,
  outf.close()  
  print "done" 
    

csvHeader='"Seq","Source","BreakSet","DatasetNo","SetNo","SetNoInData","Type","Model","Treat","%BreakSet","Year","Predecessor","Successor","%SegmentFound","Ti","Shift","%YearOverall","Modal Year","ModalValue","SecondModalYear","SecondModalValue","SegmentTrend","PreTrend",PostTrend","ShiftFromTrends","30YearTrend","15YrTi","15YearShift","30Yr Mode","30YrModePct","30YrMode2","30YrMode2%","Stability"'  
ancovaHeader=',":","Byear","slope1","Pr(slope1)","Sig(slope1)","Rsq(slope1)","slope2","Pr(slope2)","Sig(slope2)","Rsq(slope2)","non-trend-shift","Pr(non-trend-shift)","Sig(non-trend-shift)",":","PrTrends","sig(Trends)","Pr(shifts)","sig(shifts)",Pr(Regime)","sig(Regime)"'


def newstyle(csvrepfn,outfn): #from csv files
  outf=open(outfn, "w")
  with open(csvrepfn,'r') as f:
    lines=f.readlines()
  for line in lines:
    #print line
    if line[:5] == '"Seq"':
      print >>outf,csvHeader+',":",'+cmip5_info.header+ancovaHeader
    else:
      if line.find("no breaks")>=0:
        nobreaks = True
      else:
        nobreaks = False
      words=line.split('"')
      
      try:
  #      nobreaks=False
        fn=os.path.basename(words[0].split(",")[1])
      except Exception as e:
        print "Oops", line, words, str(e)
        if line.find("no breaks")>=0:
          nobreaks = True
          pass
        else:
          raise e
      #print fn+":"
      if len(words) < 2:
        print "oops2" , line, words
      if not nobreaks:
        yset=eval(words[1])
        #print "YSET", yset
        yrss=words[2].split(',')[8:11]
        yr=str(yrss[0])
        pre=str(yrss[1])
        post=str(yrss[2])
        
       # print yr,pre,post
        dfn=fn[:-len('final.csv')-1]
        if dfn[:len("Global-temps")]=="Global-temps":
          dfn = "Global temps"+dfn[len("Global-temps"):]
        tfn=os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\"+dfn
        #print fn, tfn 
        if dfn[:len("Global")]=="Global" or dfn[0] == 'r':
          print "do",tfn
          try:
            direc=glob.glob(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\"+dfn+'/*.trace')[0]
            bp=statbreaks.brkrpt(direc)
            #print bp.years()
            yrix=int(list(bp.years()).index(float(yr)))
            preix=int(list(bp.years()).index(float(pre)))
            postix=int(list(bp.years()).index(float(post)))
            #print bp.years()[preix:postix+1]
            ancova=doancova.ANCOVA(bp.years()[preix:postix+1], bp.ys()[preix:postix+1], bp.years()[yrix],Rpath='"C:\\Program Files\\R\\R-3.0.2\\bin\\Rscript.exe"')
            ancova.parseResult(ancova.results())
            print >>outf,line[:-1],',":"',
            if dfn[:len("Global")]=="Global":
              info=cmip3_info.matchfn(tfn)
              for inf in info:
                print >>outf, ",",inf,
            elif dfn[0] == 'r':
              info=cmip5_info.matchfn(tfn)
              for inf in info:
                print >>outf, ",",inf,
            print >>outf, ',":"',
            for inf in [
                ancova._ANCOVA__breakYear,
                ancova.firstSlope(),
                formatprobs(ancova.firstSlopeProb()), 
                ancova.firstRsq(),
                ancova.secondSlope(),
                formatprobs(ancova.secondSlopeProb()),
                ancova.secondRsq(),
                ancova.shift(),
                formatprobs(ancova._ANCOVA__prDiffShift),
                ":",
                formatprobs(ancova.changeOfTrendProb(), ancova._ANCOVA__prDiffShift),
                formatprobs(ancova.changeOfTrendProb()* ancova._ANCOVA__prDiffShift)
                ]:
              print >>outf,",",inf,
            print >>outf
            #if 
          except Exception as e:
            fl=os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\"+fn[:-len('final.csv')-1]+'/*.trace'
            print >>outf,line[:-1],"Error getting files",fl
            print "using ", fl
            #print glob.glob(fl)
            #raise Exception(str(e))
            pass
            
          #print direc,
  outf.close()  
  print "done" 
    
if __name__ == "__main__":
  #test case newstyle(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\cmip5rcp26_027.csv", "test.txt")    
  #newstyle(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\pitest.csv", "pitest.csv.txt")      
  #NOAA as a unspecifuied bname
  #newstyle(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\pitestNOAA19.csv", "pitest.csv.txt")      
  '''
  newstyle(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\cmip5_piControl.repeat.csv", "piControlr.txt")    
  newstyle(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\cmip5_piControl.csv", "piControl.txt")  
  
  newstyle(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\GlobaTemps.repeat.csv","GlobaTemprs.txt")
     
  newstyle(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\GlobaTemps.csv","GlobaTempr.txt")      
  newstyle(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\cmip5_RCP26.csv","RCP26.txt")      
  newstyle(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\cmip5_RCP45.csv","RCP45.txt") 
  newstyle(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\cmip5_RCP60.csv","RCP60.txt")    
  newstyle(os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\cmip5_RCP85.csv","RCP85.txt")  
  '''       
  import uniqsets
  uniqsets.doit("piControlr.txt","piControlr.txt.csv")
  uniqsets.doit("piControl.txt","piControl.txt.csv")
  uniqsets.doit("GlobaTemprs.txt","GlobaTemps.repeat.txt.csv")
  uniqsets.doit("GlobaTempr.txt","GlobaTemps.txt.csv")
  
  uniqsets.doit("RCP26.txt","RCP26.txt.csv")
  uniqsets.doit("RCP45.txt","RCP45.txt.csv")
  uniqsets.doit("RCP60.txt","RCP60.txt.csv")
  uniqsets.doit("RCP85.txt","RCP85.txt.csv")
  print "done"