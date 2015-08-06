# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 12:07:38 2015

@author: s4493222
"""
import numpy as np
import csvfile
import statbreaks
import glob
#CollateRpt
SVNRevision="$Revision: 346 $"

#input produces as
#IFS=$'\n';for f in *.final.csv; do while read line; do echo "$f $line";done <$f;done >x.x


class CollateRpt(object):
  def __init__(self, rootname):
    self.__rootname=rootname
    self.__finalfile=rootname+'.final.csv'
    self.__lines=[]
    with open(self.__finalfile) as f:
      lines=f.readlines()
      for line in lines:
        yrset=line[:line.find('"',1)+1]
        tail=line[line.find('"',1)+2:-1]
        self.__lines.append([yrset, tail.split(' ')])
      #print self.__lines[0][1][0] will be the frequency count
    #now we have a list of lists [[str(breaks), [str(parameters)]]]
    #find the first line with the max freq
    maxi=0
    maxix=0
    maxs=[]
    self.__METADATA=[]
    for l in range(len(self.__lines)):
      try:
        fr=int(self.__lines[l][1][0])
        if fr > maxi:
          maxi=fr
          maxix=l
          maxs=self.__lines[l][0]
      except:  #metadate at bottom
        self.__METADATA.append([rootname, self.__lines[l]])
        pass
    #open any trace file in the _self.__rootnamedir and extract the years and yi values
    samplefn = glob.glob(rootname+'\\*.trace')[0]
    sb=statbreaks.brkrpt(samplefn)
    self.__years=sb.years()
    self.__ys=sb.ys()
    self.__maxset=maxs
    #now finds all the lines nelonging to modal
    self.__analysis = []
    for line in self.__lines:
      if maxs == line[0]:
        self.__analysis.append(line[1])
    #parameters of interest are Year, Ti, Shift, Freq, 
    #                 line[1][1, 5, 6,0]
    
    
    #print self.__analysis
    
  def params(self):
    return [line[[1,5,6,0]] for line in np.array(self.__analysis)]
    #return [line[[1,5,6,0,-8]] for line in np.array(self.__analysis)] # gives trend to trend shift
    
  def rootname(self):
    return self.__rootname
  def finalfile(self):
    return self.__finalfile
  def lines(self):
    return self.__lines
  def years(self):
    return self.__years
  def ys(self):
    return self.__ys
  def modalbreaks(self):
    return self.__maxset
  def yearset(self):
    return set([int(y) for y in self.__years])
  def metadata(self):
    return self.__METADATA

remonth={
'Jan':'01Jan',
'Feb':'02Feb',
'Mar':'03Mar',
'Apr':'04Apr',
'May':'05May',
'Jun':'06Jun',
'Jul':'07Jul',
'Aug':'08Aug',
'Sep':'09Sep',
'Oct':'10Oct',
'Nov':'11Nov',
'Dec':'12Dec',
}

def MonthSortFileList(filesin):
  #step1, replace months with dummy names
  #sort
  #step 3 replace dummy names with names
  import copy
  files=copy.copy(filesin)
  for i in range(len( files)):
    for k in remonth:
      files[i]=files[i].replace(k, remonth[k])
  files.sort()
  for i in range(len( files)):
    for k in remonth:
      files[i]=files[i].replace(remonth[k],k)
  return files
  
  
def printCollation(filelist):
  individuals={}
  yearset = set()
  keylist=[]
  print "years,",
  for fi in MonthSortFileList(filelist):
    print ",",'"'+fi+'"',",",'"'+fi+'_Ti"',",",'"'+fi+'_shift"',",",'"'+fi+'_freq"',
    cr=CollateRpt(fi)
    keylist.append(fi)
    yearset.update(cr.yearset())
    individuals[fi]=cr
  print
  years = list(yearset)
  years.sort()
  for y in years:
    print y,",",
    for fi in keylist: #filelist 
      cr=individuals[fi]
      params=cr.params()
      if y in cr.years():
        print ","+str(cr.ys()[np.where(cr.years()==y)][0]),
        empt=True
        for p in params:
          if int(p[0]) == y:
            print ","+p[1]+","+p[2]+",",p[3],
            empt=False
        if empt:
          print ",0,0,0",
      else:
        print ",,,,",
    print    
    
  print '"years/ys",',
#  for fi in filelist:
  for fi in keylist:
    print ",",'"'+fi+'"',
  print
  for y in years:
    print y,",",
    for fi in keylist:#filelist:
      cr=individuals[fi]
      params=cr.params()
      if y in cr.years():
        print ","+str(cr.ys()[np.where(cr.years()==y)][0]),
      else:
        print ",",
    print    

  print '"years/Ti",',
#  for fi in filelist:
  for fi in keylist:
    print ",",'"'+fi+'_Ti"',
  print
  for y in years:
    print y,",",
    for fi in keylist:
      cr=individuals[fi]
      params=cr.params()
      empt=True
      if y in cr.years():
        for p in params:
          if int(p[0]) == y:
            print ","+p[1],
            empt=False
        if empt:
          print ",0",
      else:
        print ", ",
    print    

  print '"years/Shifts",',
#  for fi in filelist:
  for fi in keylist:
    print ",",'"'+fi+'_Ti"',
  print
#  years = list(yearset)
#  years.sort()
  for y in years:
    print y,",",
    for fi in keylist:#filelist:
      cr=individuals[fi]
      params=cr.params()
      empt=True
      if y in cr.years():
        for p in params:
          if int(p[0]) == y:
            print ","+p[2],
            empt=False
        if empt:
          print ",0",
      else:
        print ", ",
    print    

  print '"years","freq",'
#  for fi in filelist:
  for fi in keylist:
    cr=individuals[fi]
    params=cr.params()
    try:
      print fi, ",,",params[0][3]
    except:
      print ",N/A"
      pass
    
  print "metadata"
#  for fi in filelist:
  for fi in keylist:
    cr=individuals[fi]
    print fi, 
    md=cr.metadata()
    for m in md:
      try:
        print m
      except:
        print m
'''
adjusted.csv_cru.final.csv "[1979, 1987, 2010]" 5 1987 1979 2010 5 24.8326335798 0.292969542243 99 1996 0.91 1994 0.08 0.0170225424609 0.0134845833333 0.0174906702899 0.00600804347842 0.0145624142157 9.25247370132 0.121960599766 1987 0.79 1986 0.12
adjusted.csv_cru.final.csv "[1979, 1986, 2010]" 1 1986 1979 2010 1 24.803505973 0.292355192605 1 1996 0.95 1994 0.03 0.0170225424609 0.0140699404762 0.017605384058 -0.00264160628018 0.0134227573529 9.87283521925 0.134764915001 1982 0.73 1987 0.19
adjusted.csv_cru.final.csv "[1979, 1987, 1996, 2010]" 94 1987 1979 1996 94 11.1896300458 0.132740979242 99 1987 0.84 1986 0.1 0.0139769865841 0.0134845833333 0.0103895833333 0.0259507870373 0.0145624142157 9.33699477961 0.121901958535 1987 0.91 1986 0.07
adjusted.csv_cru.final.csv "[1979, 1987, 1996, 2010]" 94 1996 1987 2010 94 18.9333459287 0.234727189531 94 1996 0.94 1997 0.05 0.017605384058 0.0121030808081 0.0122833150183 0.0884574642025 0.0207895343137 12.1595602276 0.186448433688 1996 0.99 1993 0.01
adjusted.csv_cru.final.csv closed  18:59:11,Wed,18-Feb-2015
adjusted.csv_cru.final.csv Version numbers
adjusted.csv_cru.final.csv shuffleRpt $Revision: 346 $
adjusted.csv_cru.final.csv shuffle $Revision: 346 $
adjusted.csv_cru.final.csv bivariate_multi as bivariate $Revision: 346 $
adjusted.csv_cru.final.csv ConvergentBreaks as convergent_breaks $Revision: 346 $
adjusted.csv_cru.final.csv regress $Revision: 346 $
adjusted.csv_cru.final.csv statbreaks $Revision: 346 $
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv "[1861.0, 1927, 1964, 1977, 1995, 2015, 2029, 2079, 2100.0]" 7 1927 1861.0 1964 100 33.599049032 0.0918385799604 100 1927 0.9 1943 0.08 0.00112755507665 -0.000775170272298 0.00265517622741 0.0670325697361 0.00788415439756 5.97937117336 0.126148367841 1922 0.36 1923 0.3
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv "[1861.0, 1927, 1964, 1977, 1995, 2015, 2029, 2079, 2100.0]" 7 1964 1927 1977 100 28.4041262691 0.169137259408 100 1964 1.0 0 0 0.00505985294067 0.00278323667614 0.00750636444793 0.0653065285644 0.00958702206502 9.2493047153 0.0984876179478 1964 0.95 1965 0.05
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv "[1861.0, 1927, 1964, 1977, 1995, 2015, 2029, 2079, 2100.0]" 7 1977 1964 1995 100 21.0452795929 0.20084659057 100 1977 0.84 1978 0.1 0.0116113712752 0.0103534798467 0.00832792398594 0.0542658729952 0.0143193749875 8.90614147619 0.12976686176 1977 0.86 1975 0.08
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv "[1861.0, 1927, 1964, 1977, 1995, 2015, 2029, 2079, 2100.0]" 7 1995 1977 2015 7 29.0760169625 0.320612067687 100 1995 1.0 0 0 0.014578689269 0.00992270467849 0.00800325187749 0.147927938601 0.0215314583195 9.90612324993 0.20510296877 1995 0.49 2000 0.28
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv "[1861.0, 1927, 1964, 1977, 1995, 2015, 2029, 2079, 2100.0]" 7 2015 1995 2029 7 19.5198115605 0.255026471945 7 2015 0.63 2012 0.15 0.0140408986919 0.0108818722991 0.0187149450066 0.00584855457828 0.0208938970778 9.45867593644 0.20325361602 2011 0.78 2012 0.16
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv "[1861.0, 1927, 1964, 1977, 1995, 2015, 2029, 2079, 2100.0]" 7 2029 2015 2079 7 45.6602745681 0.322332958731 100 2029 0.94 2033 0.04 0.00662681927214 0.0211513690189 0.00237622168496 0.113898319434 0.0161203921882 8.45486225041 0.142138749056 2029 0.88 2025 0.08
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv "[1861.0, 1927, 1964, 1977, 1995, 2015, 2029, 2079, 2100.0]" 7 2079 2029 2100.0 100 31.1754653751 -0.156233139363 100 2079 1.0 0 0 -0.00176211626526 0.0028956538433 -0.00418388528557 -0.182435227773 -0.0118550122442 10.7976707751 -0.139732110386 2079 0.98 2078 0.02
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv "[1861.0, 1927, 1964, 1977, 1995, 2019, 2029, 2079, 2100.0]" 93 1927 1861.0 1964 100 33.8356913228 0.0923555121054 100 1927 0.89 1943 0.09 0.00112755507665 -0.000775170272298 0.00265517622741 0.0670325697361 0.00788415439756 5.92334226058 0.11477258301 1923 0.51 1922 0.24
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv "[1861.0, 1927, 1964, 1977, 1995, 2019, 2029, 2079, 2100.0]" 93 1964 1927 1977 100 28.4019358536 0.169575979919 100 1964 0.98 1953 0.01 0.00505985294067 0.00278323667614 0.00750636444793 0.0653065285644 0.00958702206502 9.26671276817 0.0989041026471 1964 0.93 1965 0.05
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv "[1861.0, 1927, 1964, 1977, 1995, 2019, 2029, 2079, 2100.0]" 93 1977 1964 1995 100 21.0603166031 0.201288813633 100 1977 0.9 1978 0.07 0.0116113712752 0.0103534798467 0.00832792398594 0.0542658729952 0.0143193749875 9.02993697246 0.130835776386 1977 0.78 1975 0.12
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv "[1861.0, 1927, 1964, 1977, 1995, 2019, 2029, 2079, 2100.0]" 93 1995 1977 2019 93 32.3936028396 0.343915069075 100 1995 0.92 1999 0.06 0.0142982658826 0.00992270467849 0.00924588768163 0.138957167878 0.0215314583195 9.88169447874 0.206066413607 1995 0.41 2000 0.36
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv "[1861.0, 1927, 1964, 1977, 1995, 2019, 2029, 2079, 2100.0]" 93 2019 1995 2029 93 19.6701192969 0.256970377893 93 2015 0.66 2012 0.11 0.0140408986919 0.0112091153897 0.0172342423417 0.0437609490643 0.0209349877311 9.86359654277 0.193774807887 2021 0.89 2025 0.07
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv "[1861.0, 1927, 1964, 1977, 1995, 2019, 2029, 2079, 2100.0]" 93 2029 2019 2079 93 39.4494073803 0.252295052593 100 2033 0.8 2029 0.09 0.00550312092181 0.02067537873 0.00237622168496 0.11620151402 0.0161203921882 8.57905755731 0.142320558016 2029 0.86 2025 0.09
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv "[1861.0, 1927, 1964, 1977, 1995, 2019, 2029, 2079, 2100.0]" 93 2079 2029 2100.0 100 30.9874470883 -0.155785894783 100 2079 1.0 0 0 -0.00176211626526 0.0028956538433 -0.00418388528557 -0.182435227773 -0.0118550122442 10.9675289734 -0.141549265604 2079 0.97 2078 0.03
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv closed  19:41:35,Fri,20-Feb-2015
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv Version numbers
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv shuffleRpt $Revision: 346 $
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv shuffle $Revision: 346 $
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv bivariate_multi as bivariate $Revision: 346 $
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv ConvergentBreaks as convergent_breaks $Revision: 346 $
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv regress $Revision: 346 $
r1i1p1_FGOALS_g2_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_030.dat.final.csv statbreaks $Revision: 346 $
'''

def printFinal_icmip5(filelist, outname="output.txt"):
#input produced as
#for f in *.final.csv; do while read line; do echo "$f $line";done <$f;done >x.x
#grep r1i1p1 x.x >x.y
  header='"Seq","Source","BreakSet","DatasetNo","SetNo","SetNoInData","Repeat","Model","RCP","%BreakSet","Year","Predecessor","Successor","%Segment","Found","Ti,Shift","%YearOverall","Modal Year","ModalValue","SecondModalYear","SecondModalValue","SegmentTrend","PreTrend",PostTrend","ShiftFromTrends","15YearTrend","15YrTi","15YearShift","15Yr Mode","15YrModePct","15YrMode2","15YrMode2%"'

  for f in filelist:
    with open(f,'r') as fi:
      lines=fi.readlines()
  metadata=[]
  outf = open(outname,"w") 
  seq=0
  datasetno=1
  setno=1
  setnoindata=1
  oldfilename = None
  oldsetstr=None
  EOF=False
  BrkCount = 0
  files={}
  for line in lines:
#    if EOF and BrkCount == 0: #if no lines printed the 
#      setnoindata=0
#      print >>outf,header
#      oldfilename = filename
      
    line = line[:-1]
    lfn = line.split()[0]
    if not lfn in files:
      files[lfn]=0
    if line.find("closed") >= 0:
      metadata.append(line)
      EOF=True
      #print line
    elif (line.find("Version") >= 0):
      metadata.append(line)
    elif (line.find("$Rev") >= 0):
      metadata.append(line)
    else:
      EOF=False
      BrkCount+=1
      seq +=1
      filename = line[:line.find('"')-1]
      if oldfilename == None:
        print >>outf,header
      if oldfilename != None and oldfilename != filename:
        datasetno +=1
        setnoindata=0
        print >>outf,header
      oldfilename = filename
      files[filename]+=1
      pos1=line.find('"',1)
      pos2=line.find('"',pos1 + 1)
      setstr= line[pos1:pos2+1]
      if oldsetstr != None and oldsetstr != setstr:
        setno +=1
        setnoindata +=1
      oldsetstr = setstr
      
      tail = line[pos2+2:]
      us1=filename.find("_",1)
      rep=filename[:us1]
      rcp1=filename.find("RCP")
      if rcp1 == -1:
        rcp1=filename.find("AR5")
      icmppos=filename.find("_icmip")
      rcp=filename[rcp1:icmppos]
      model=filename[us1+1:rcp1-1]
      outline=str(seq)+","+filename+','+setstr+','+str(datasetno)+','+str(setno)+','+str(setnoindata)+','+rep+','+model+','+rcp
      for s in tail.split(' '):
        outline += ','+s
      print >>outf,outline
  for fn in files.keys():
    if files[fn] == 0:
      print >>outf, fn, " no breaks"
      print fn, " no breaks"
  print outname, "Done"
      
def parsefilename(filename):
  if filename[:len('GISSTEMP')]=='GISSTEMP':
    model="GISSTEMPto-2013b"
    rep="Observed"
    rcp=filename[len(model)+2:len(model)+8]
  elif filename[:len('adjusted.csv')] == 'adjusted.csv':
    model=filename.split("_")[1]
    if model.find('res') >-1:
      rep='F&R'
      model=model[:-3]
      rcp='res'
    elif model.find('raw') >-1:
      rep='F&R'
      model=model[:-3]
      rcp='raw'
    elif model.find('fit') >-1:
      rep='F&R'
      model=model[:-3]
      rcp='fit'
    elif model=='meanmodel':
      rep='F&R'
      model='meanmodel'
      rcp='mean'
    else:
      rep='F&R'
      model=model
      rcp="result"
  elif filename[:len('aravg.mon.land_ocean')] == 'aravg.mon.land_ocean':
    model='aravg.mon.land_ocean'
    rep="Observed"
    rcp=filename[len(model)+1:len(model)+8]
    model="NCDC-banded"
  elif filename[:len('had4_krig_annual_v2_0_0')] == 'had4_krig_annual_v2_0_0':
    model='had4_krig_annual_v2_0_0'
    rep="Observed"
    rcp='COWTAN_WAY'
  elif filename[:len('HadCRUT.4.3')] == 'HadCRUT.4.3':
    model='HadCRUT.4.3'
    rep="Observed"
    rcp=filename[filename.find('_'):]
  elif filename[:len('HadCRUT.4.2')] == 'HadCRUT.4.2':
    model='HadCRUT.4.2'
    rep="Observed"
    rcp=filename[filename.find('_'):]
  elif filename[:len('Land_and_Ocean_summary')] == 'Land_and_Ocean_summary':
    model='BERKLEY'
    rep="Observed"
    rcp="BEST"
  elif filename[:len('Arty')] == 'Arty':
    model = "TEST"
    rep = "synthetic"
    rcp = filename[-13:-9]
  elif filename[0] == "r":
    model=filename[filename.find("_"):filename.find("_",2)+1]
    rep=filename[:filename.find("_")+1]
    rcp=filename[filename.find("RCP"):filename.find("icmip")+1]
  elif filename[:3] == "tas":
    model=filename[filename.find("_"):filename.find("_",2)+1]
    rep=filename[:filename.find("_")+1]
    rcp=filename[filename.find("RCP"):filename.find("icmip")+1]
  elif filename[:len("RSS_TS")] == "RSS_TS":
    rcp="LandOcean"
    words=filename.split("_")
    model="RSS_"+words[3]
    rep=words[4]
  #"uahncdc_tp_6.0beta1_SoPolOcean.final.csv"
  elif filename[:len("uahncdc_")]=="uahncdc_":
    words=filename.split("_")
    model="UAH_"+words[1]
    words=words[-1].split(".")
    if words[1][:-len("Ocean")] == "Ocean":
      rcp="Ocean"
    elif words[1][:-len("Land")] == "Land":
      rcp="Land"
    else:
      rcp="LandOcean"
    rep = words[1]

  else: raise Exception ("Cannot process ")   
  return model, rep, rcp      
      
def printFinal_current(filelist, outname, cmip5=False):
#input produced as
#for f in *.final.csv; do while read line; do echo "$f $line";done <$f;done >x.x
#grep r1i1p1 x.x >x.y
  header='"Seq","Source","BreakSet","DatasetNo","SetNo","SetNoInData","Type","Model","Treat","%BreakSet","Year","Predecessor","Successor","%SegmentFound","Ti","Shift","%YearOverall","Modal Year","ModalValue","SecondModalYear","SecondModalValue","SegmentTrend","PreTrend",PostTrend","ShiftFromTrends","30YearTrend","15YrTi","15YearShift","30Yr Mode","30YrModePct","30YrMode2","30YrMode2%","Stability"'
  header='"Seq","Source","BreakSet","DatasetNo","SetNo","SetNoInData","Type","Model","Treat","%BreakSet","Year","Predecessor","Successor","%SegmentFound","Ti","Shift","%YearOverall","ModalYear","ModalValue","SecondModalYear","SecondModalValue","SegmentTrend","PreTrend",PostTrend","ShiftFromTrends","30YearTrend","30YrTi","30YearShift","30Yr Mode","30YrModevalue","30YrMode2","30YrMode2Value","Stability"'
  if cmip5:
    header='"Seq","Source","BreakSet","DatasetNo","SetNo","SetNoInData","Repeat","Model","RCP","%BreakSet","Year","Predecessor","Successor","%SegmentFound","Ti","Shift","%YearOverall","ModalYear","ModalValue","SecondModalYear","SecondModalValue","SegmentTrend","PreTrend",PostTrend","ShiftFromTrends","30YearTrend","30YrTi","30YearShift","30Yr Mode","30YrModevalue","30YrMode2","30YrMode2Value","Stability"'

  for f in filelist:
    with open(f,'r') as fi:
      lines=fi.readlines()
  metadata=[]
  outf = open(outname,"w") 
  seq=0
  datasetno=1
  setno=1
  setnoindata=1
  oldfilename = None
  oldsetstr=None
  files={}

  for line in lines:
    subdone=False
    line = line[:-1]
    if line[:len('Global temps')]=='Global temps':
      line=line[:6]+"-"+line[7:]
      subdone=True
    lfn = line.split()[0]
    if not lfn in files:
      files[lfn]=0
    if line.find("closed") >= 0:
      metadata.append(line)
    elif (line.find("Version") >= 0):
      metadata.append(line)
    elif (line.find("$Rev") >= 0):
      metadata.append(line)
    else:
      seq +=1
      filename = line[:line.find('"')-1]
      if subdone:
        filename=filename[:6]+"-"+filename[7:]
      
      files[filename]+=1
      if oldfilename == None:
        print >>outf,header
    
      if oldfilename != None and oldfilename != filename:
        datasetno +=1
        setnoindata=0
        print >>outf,header
      #print oldfilename, filename
      oldfilename = filename
      
      if filename[:len('GISSTEMP')]=='GISSTEMP':
        model="GISSTEMPto-2013b"
        rep="Observed"
        rcp=filename[len(model)+2:len(model)+8]
      elif filename[:len('adjusted.csv')] == 'adjusted.csv':
        model=filename.split("_")[1]
        if model.find('res') >-1:
          rep='F&R'
          model=model[:-3]
          rcp='res'
        elif model.find('raw') >-1:
          rep='F&R'
          model=model[:-3]
          rcp='raw'
        elif model.find('fit') >-1:
          rep='F&R'
          model=model[:-3]
          rcp='fit'
        elif model=='meanmodel':
          rep='F&R'
          model='meanmodel'
          rcp='mean'
        else:
          rep='F&R'
          model=model
          rcp="result"
      elif filename[:len('aravg.mon.land_ocean')] == 'aravg.mon.land_ocean':
        model='aravg.mon.land_ocean'
        rep="Observed"
        rcp=filename[len(model)+1:len(model)+8]
        model="NCDC-banded"
      elif filename[:len('had4_krig_annual_v2_0_0')] == 'had4_krig_annual_v2_0_0':
        model='had4_krig_annual_v2_0_0'
        rep="Observed"
        rcp='COWTAN_WAY'
      elif filename[:len('HadCRUT.4.3')] == 'HadCRUT.4.3':
        model='HadCRUT.4.3'
        rep="Observed"
        rcp=filename[filename.find('_'):]
      elif filename[:len('HadCRUT.4.2')] == 'HadCRUT.4.2':
        model='HadCRUT.4.2'
        rep="Observed"
        rcp=filename[filename.find('_'):]
      elif filename[:len('Land_and_Ocean_summary')] == 'Land_and_Ocean_summary':
        model='BERKLEY'
        rep="Observed"
        rcp="BEST"
      elif filename[:len('Arty')] == 'Arty':
        model = "TEST"
        rep = "synthetic"
        rcp = filename[-13:-9]
#      elif filename[0] == "r":
#        model=filename[filename.find("_"):filename.find("_",2)+1]
#        rep=filename[:filename.find("_")+1]
#        rcp=filename[filename.find("RCP"):filename.find("icmip")+1]
      elif filename[:3] == "tas":
        model=filename[filename.find("_"):filename.find("_",2)+1]
        rep=filename[:filename.find("_")+1]
        rcp=filename[filename.find("RCP"):filename.find("icmip")+1]
#chnages for the Hadley HADCRUT annual and monthl plus CRUT plus HadSST which came via CRU
      elif filename[:len('CRUTEM4-')] == 'CRUTEM4-':
        model='CRUTEM4'
        posn1=filename.find("_")
        if posn1 <0:
          rcp="Annual"
        else:
          rcp=filename[posn1+1:posn1+4]
        rep=filename[len('CRUTEM4-'):len('CRUTEM4-')+2]
      elif filename[:len('HadCRUT4-')] == 'HadCRUT4-':
        model='HadCRUT4'
        posn1=filename.find("_")
        if posn1 <0:
          rcp="Annual"
        else:
          rcp=filename[posn1+1:posn1+4]
        rep=filename[len('HadCRUT4-'):len('HadCRUT4-')+2]
      elif filename[:len('HadSST3-')] == 'HadSST3-':
        model='HadSST3'
        posn1=filename.find("_")
        if posn1 <0:
          rcp="Annual"
        else:
          rcp=filename[posn1+1:posn1+4]
        rep=filename[len('HadSST3-'):len('HadSST3-')+2]
#chnages NCDC Monthly set, Llan, Ocean, LandOcean
      elif filename.find("_LandOcean-") >=0:
        model='NCDCLandOcean'
        posn1 = filename.find("_LandOcean-")+len("_LandOcean-")-1
        if filename.find("Jan-Dec") >0:
          rcp="Annual"
        else:
          rcp=filename[posn1+1:posn1+4]
        rep=filename[:2]
      elif filename.find("-Land-") >=0:
        model='NCDCLand'
        posn1 = filename.find("-Land-")+len("-Land-")-1
        if filename.find("Jan-Dec") >0:
          rcp="Annual"
        else:
          rcp=filename[posn1+1:posn1+4]
        rep=filename[:2]
      elif filename.find("-Ocean-") >=0:
        model='NCDCOcean'
        posn1 = filename.find("-Ocean-")+len("-Ocean-")-1
        if filename.find("Jan-Dec") >0:
          rcp="Annual"
        else:
          rcp=filename[posn1+1:posn1+4]
        rep=filename[:2]
#The HadSST.3.1 monthly analysis        
      elif filename.find("HadSST.3.1") >=0:
        model='HadSST.3.1.1.0'
        posn1 = filename.find("_ts_")+len("_ts_")+1
        if filename.find("ts.final") >0:
          rcp="Annual"
          posn1=filename.find("_ts.")+len("_ts.")+1
        else:
          rcp=filename[posn1-1:posn1+2]
        posn2=filename.find("_monthly_")+len("_monthly_")+1
        rep=filename[posn2-1:posn1-5]
        
#the cmip5 resiulst downloaded with names commencing with r1i1p1 etc
      elif filename[:filename.find("_")].translate(None,"0123456789") == "rip":
        if filename.find("piControl") >= 0:
          rcp = "piControl"
          words= filename.split("_")
          rep=words[0]
          model=words[1]
          if model=="FGOALS" and words[2] == "g2":
            model="FGOALS_g2"#x.x_6May
        else:
          words= filename.split("_")
          rep=words[0]
          model=words[1]
          rcp=words[2]
          if model=="FGOALS" and rcp == "g2":
            model="FGOALS_g2"
            rcp=words[3]
      elif filename[:len('Global-temps')]=='Global-temps':
        tail=filename[filename.find('_')+1:filename.find('.',15)]
        if tail[-1]=="E" : #@ensemble - no rep
          rep="Ens"
          tail=tail[:-1]
        else:
          rep=tail[tail.find('r')+1:]
          tail=tail[:tail.find('r')]
        pos=tail.find('A1B')
        if pos >= 0:
          rcp='SRES_A1B'
        else:
          pos=tail.find('A1')
          if pos >= 0:
            rcp='SRES_A1'
          else:
            pos=tail.find('A2')
            if pos < 0:
              raise "Global-temps what SRES", tail
            else:
              rcp='SRES_A2'
        model=tail[:pos]
      #RSS and UAH JHR 6May
      elif filename[:len("RSS_TS")] == "RSS_TS":
        rcp="LandOcean"
        words=filename.split("_")
        model="RSS_"+words[3]
        rep=words[4]
      #"uahncdc_tp_6.0beta1_SoPolOcean.final.csv"
      elif filename[:len("uahncdc_")]=="uahncdc_":
        words=filename.split("_")
        model="UAH_"+words[1]
        words=words[-1].split(".")
        if words[1][:-len("Ocean")] == "Ocean":
          rcp="Ocean"
        elif words[1][:-len("Land")] == "Land":
          rcp="Land"
        else:
          rcp="LandOcean"
        rep = words[1]
      else: raise Exception ("Cannot process "+filename)   
      
      pos1=line.find('"',1)
      pos2=line.find('"',pos1 + 1)
      setstr= line[pos1:pos2+1]
      if oldsetstr != None and oldsetstr != setstr:
        setno +=1
        setnoindata +=1
      oldsetstr = setstr
      tail = line[pos2+2:]
      

      outline=str(seq)+","+filename+','+setstr+','+str(datasetno)+','+str(setno)+','+str(setnoindata)+','+rep+','+model+','+rcp
      for s in tail.split(' '):
        outline += ','+s
      print >>outf,outline
  for fn in files.keys():
    if files[fn] == 0:
      print >>outf, fn, " no breaks"
      print fn, " no breaks"
      model, rep, rcp = parsefilename(filename)
      outline="no breaks"+","+filename+','"N/A"+','+"N/A"+','+"N/A"+','+"N/A"+','+rep+','+model+','+rcp+",,,,,,,"
      print >>outf, outline
  print outname, "Done"
      
      
    
#  
if __name__ == "__main__":
#  printFinal_current(['NOAA.x'], "NOAA.csv",cmip5=True)
#  printFinal_current(['cmip5_piControl.x'], "cmip5_piControl.csv",cmip5=True)
#  printCollation(glob.glob("r*i*p*piControl*[0123456789]"))  
#  printFinal_current(['cmip5_RCP26.x'], "cmip5_RCP26.csv",cmip5=True)
#  printCollation(glob.glob("r*i*p*_*_RCP*[0123456789].dat"))  
#  printFinal_current(['cmip5_RCP45.x'], "cmip5_RCP45.csv",cmip5=True)
#  printCollation(glob.glob("r*i*p*_*_*RCP4.5_icmip*[0123456789]"))  
#  printFinal_current(['cmip5_RCP60.x'], "cmip5_RCP60.csv",cmip5=True)
#  printCollation(glob.glob("r*i*p*_*_*RCP6_icmip*[0123456789]"))  
#  printFinal_current(['cmip5_RCP85.x'], "cmip5_RCP85.csv",cmip5=True)
#  printCollation(glob.glob("r*i*p*_*_*RCP8.5_icmip*[0123456789]"))  
#  printFinal_current(['GlobaTemps.x'], "GlobaTemps.csv",cmip5=True)
#  printCollation(glob.glob("Global*temps*[0123456789E]"))  
#  printFinal_current(['GISS15Apr.x'], "GISS15Apr.csv",cmip5=True)
#GISS15Apr.x
#  printCollation(glob.glob("GISSTEMPv3-April2015*[NUSmb]"))  
#  printFinal_current(['GISSDIFFS.x'], "GISSDIFFS.csv",cmip5=True)
#GISS15Apr.x
#  printCollation(glob.glob("GISSTEMP_v3_April_2015_withDiffs_from_2014*[NUSmb]"))  

#6May   
  printFinal_current(['rss.x'], "RSS.csv",cmip5=False)
  printCollation(glob.glob("RSS_TS*[TS]")) 
  
#  printFinal_current(['uah.x'], "UAH.csv",cmip5=False)
#  printCollation(glob.glob("uahncdc_*[lnd8tH9T]")) 
  '''
  printFinal_current(['GISSTEMP.x'], "GISSTEMP2.csv")
  printFinal_current(['NCDC_Banded.x'], "NCDC_Banded2.csv")
  printFinal_current(['COWTAN_WAY.x'], "COWTAN_WAY2.csv")
  printFinal_current(['HadCRUT.4.3.x'], "HadCRUT2.4.3.csv")
  printFinal_current(['HadCRUT.4.2.x'], "HadCRUT2.4.2.csv")
  printFinal_current(['BERKLEY.x'], "BERKLEY2.csv")

  printFinal_current(['x.x_9mar'], "x.x_9mar.csv")
  printFinal_current(['HADCRUT4_CRU.x'], "HADCRUT4_CRU.csv")
  printFinal_current(['HADSST3_CRU.x'], "HADSST3_CRU.csv")
  printFinal_current(['CRUTEM4_CRU.x'], "CRUTEM4_CRU.csv")
  printFinal_current(['NCDCLandOcean.x'], "NCDCLandOcean.csv")
  printFinal_current(['NCDCOcean.x'], "NCDCOcean.csv")
  printFinal_current(['NCDCLand.x'], "NCDCLand.csv")
  printFinal_current(['HadSST.3.1.1.0.x'], "HadSST.3.1.1.0.csv")
  printFinal_current(['NCDCLandOcean.x'], "NCDCLandOcean.csv")
  printFinal_current(['NCDCOcean.x'], "NCDCOcean.csv")
  printFinal_current(['NCDCLand.x'], "NCDCLand.csv")
  '''
#  printCollation(glob.glob("HadSST.3.1.1.0_monthly_*_ts"))  
#  printCollation(glob.glob("HadSST.3.1.1.0_monthly_*_ts_???"))  
#  printCollation(glob.glob("*-Land-Jan-Dec-1880-2015"))  
#  printCollation(glob.glob("*-Land-???-1880-2015"))  
#  printCollation(glob.glob("*LandOcean-Jan-Dec-1880-2015"))  
#  printCollation(glob.glob("*LandOcean-???-1880-2015"))  
#  printCollation(glob.glob("*-Ocean-Jan-Dec-1880-2015"))  
#  printCollation(glob.glob("*-Ocean-???-1880-2015"))  
#  printCollation(glob.glob("HadSST3*.dat_???"))  
#  printCollation(glob.glob("HadSST3*.dat"))  
#  printCollation(glob.glob("HadSST3*.dat_???"))  
#  printCollation(glob.glob("CRUTEM4*.dat"))  
#  printCollation(glob.glob("CRUTEM4*.dat_???"))  

#  printCollation(glob.glob("HadCRUT4*.dat"))  
#  printCollation(glob.glob("HadSST3*.dat_???"))  
#  printCollation(glob.glob("HadSST3*.dat"))  
#  printCollation(glob.glob("HadSST3*.dat_???"))  
#  printCollation(glob.glob("CRUTEM4*.dat"))  
#  printCollation(glob.glob("CRUTEM4*.dat_???"))  
  '''  
  printFinal_current(['GISSTEMP.x'], "GISSTEMP_3.csv")
  printFinal_current(['NCDC_Banded.x'], "NCDC_Banded_3.csv")
  printFinal_current(['COWTAN_WAY.x'], "COWTAN_WAY_3.csv")
  '''
  #printFinal_current(['HadCRUT4.3.x'], "HadCRUT4.3_3.csv")
  #printFinal_current(['HadCRUT4.2.x'], "HadCRUT4.2_3.csv")
  '''
  printFinal_current(['BERKLEY.x'], "BERKLEY_3.csv")
  printFinal_current(['ARTIFICIAL.x'], "ARTIFICIAL_3.csv")
  printFinal_current(['FR.x'], "FR_3.csv")

  printFinal_current(['x.14_mar'], "x.14_mar.csv")
   '''
#   printFinal_icmip5(['piControl.x'], "piControl-2.csv")
#   printFinal_icmip5(['RCP26.x'], "RCP26-2.csv")
#   printFinal_icmip5(['RCP45.x'], "RCP45-2.csv")
#   printFinal_icmip5(['RCP6.x'], "RCP6-2.csv")
#   printFinal_icmip5(['RCP85.x'], "RCP85-2.csv")

#  cr=CollateRpt("r6i1p1_CSIRO-Mk3-6-0_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_023.dat")
#  printCollation(glob.glob("GISSTEMPto*[NUS]"))
#  printCollation(glob.glob("aravg*8"))#
#  printCollation(glob.glob("had4_krig_annual_v2_0_0"))
#  printCollation(["HadCRUT.4.2.0.0.annual_30S_30N","HadCRUT.4.2.0.0.annual_nh","HadCRUT.4.2.0.0.annual_ns_avg","HadCRUT.4.2.0.0.annual_sh"])
#  printCollation(["HadCRUT.4.3.0.0.annual_30S_30N","HadCRUT.4.3.0.0.annual_nh","HadCRUT.4.3.0.0.annual_ns_avg","HadCRUT.4.3.0.0.annual_sh"])
#  printCollation(["Land_and_Ocean_summary"])
#  printCollation(glob.glob("adjusted*[utwslch]"))
#  printCollation(glob.glob("arty*[ABCD][012345]"))
  '''
  lst=glob.glob("GISSTEMPto*[NUS]")
  lst.extend(glob.glob("aravg*8"))

  lst.extend(["had4_krig_annual_v2_0_0","HadCRUT.4.2.0.0.annual_30S_30N","HadCRUT.4.2.0.0.annual_nh","HadCRUT.4.2.0.0.annual_ns_avg","HadCRUT.4.2.0.0.annual_sh","HadCRUT.4.3.0.0.annual_30S_30N","HadCRUT.4.3.0.0.annual_nh","HadCRUT.4.3.0.0.annual_ns_avg","HadCRUT.4.3.0.0.annual_sh","Land_and_Ocean_summary"])

  
  printCollation(lst)
  '''
#  printFinal_current(['GISSTEMP.x'], "GISSTEMP.csv")
#  printFinal_current(['NCDC_Banded.x'], "NCDC_Banded.csv")
#  printFinal_current(['COWTAN_WAY.x'], "COWTAN_WAY.csv")
#  printFinal_current(['HadCRUT.4.3.x'], "HadCRUT.4.3.csv")
#  printFinal_current(['HadCRUT.4.2.x'], "HadCRUT.4.2.csv")
#  printFinal_current(['BERKLEY.x'], "BERKLEY.csv")
#  printFinal_current(['FosterRams.x'], "FosterRams.csv")

 
#  rcp26='r1i1p1*RCP2.6_icmip5*.dat'
#  rcp45='r1i1p1*RCP4.5_icmip*[0123456789]'
#  rcp6='r1i1p1*RCP6*_icmip*[0123456789]'
#  rcp85='r1i1p1*RCP8*_icmip*[0123456789]'
#  printFinal_icmip5(['xz.y'])
#  printCollation(glob.glob(rcp85))
#  for fi in glob.glob(rcp26):
#    print fi, CollateRpt(fi).params()
#    
#  for fi in glob.glob(rcp45):
#    print fi, CollateRpt(fi).params()
#
#  for fi in glob.glob(rcp6):
#    print fi, CollateRpt(fi).params()
#  for fi in glob.glob(rcp85):
#    print fi, CollateRpt(fi).params()