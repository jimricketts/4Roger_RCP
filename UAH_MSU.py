# -*- coding: utf-8 -*-
"""
Created on Mon May 04 15:44:26 2015

@author: s4493222
"""

#UAH format.

UAH_MasterDir="C:\\Users\\s4493222\\Documents\\ReferenceData\UAH\V6\\"
SVNRevision="$Revision: 346 $"
def makedatafilenames():
  global UAH_TLTDir
  global UAH_TLSDir
  global UAH_TMTDir
  global UAH_TTPDir
  global dataset
  
  UAH_TLTDir=UAH_MasterDir+"tlt\\"
  UAH_TLSDir=UAH_MasterDir+"tls\\"
  UAH_TMTDir=UAH_MasterDir+"tmt\\"
  UAH_TTPDir=UAH_MasterDir+"ttp\\"
  dataset={
  "tlt":(UAH_TLTDir,"uahncdc_lt_6.0beta1.txt"),
  "tls":(UAH_TLSDir,"uahncdc_ls_6.0beta1.txt"),
  "tmt":(UAH_TMTDir,"uahncdc_mt_6.0beta1.txt"),
  "ttp":(UAH_TTPDir,"uahncdc_tp_6.0beta1")
  }

UAH_TLTDir=UAH_MasterDir+"tlt\\"
UAH_TLSDir=UAH_MasterDir+"tls\\"
UAH_TMTDir=UAH_MasterDir+"tmt\\"
UAH_TTPDir=UAH_MasterDir+"ttp\\"
dataset={
"tlt":(UAH_TLTDir,"uahncdc_lt_6.0beta1.txt"),
"tls":(UAH_TLSDir,"uahncdc_ls_6.0beta1.txt"),
"tmt":(UAH_TMTDir,"uahncdc_mt_6.0beta1.txt"),
"ttp":(UAH_TTPDir,"uahncdc_tp_6.0beta1.txt")
}

makedatafilenames()

HeaderPrototype="Year Mo Globe  Land Ocean   NH   Land Ocean   SH   Land Ocean Trpcs  Land Ocean NoExt  Land Ocean SoExt  Land Ocean NoPol  Land Ocean SoPol  Land Ocean USA48 USA49  AUST"
Headers="Year Mo Globe  GlobeLand GlobeOcean   NH   NHLand NHOcean   SH   SHLand SHOcean Trpcs  TrpcsLand TrpcsOcean NoExt  NoExtLand NoExtOcean SoExt  SoExtLand SoExtOcean NoPol  NoPolLand NoPolOcean SoPol  SoPolLand SoPolOcean USA48 USA49  AUST".split()
import numpy as np

class UAH_MSU(object):
  def __init__(self,msulevel):
    (path,basefn)=dataset[msulevel]
    self.__filename=path+basefn
    with open(self.__filename,'r') as f:
      lines=f.readlines()
    HeaderIn=False
    SecondHeader=False
    Year1=True
    self.__metadata=[]
    self.__data=[]
    self.__dates={}
    self.__lofull=0
    self.__hifull=0
    
    lineno =0
    for line in lines: #lead blank, trailing lf
      #print len(line),line[:-1]
      line = line[1:-1]
      if len(line) > 0:
        if line ==  HeaderPrototype:
          if not HeaderIn:
            HeaderIn = True
          else:
            SecondHeader = True
        elif SecondHeader:
          self.__metadata.append(line)
        else:
          words=line.replace("-"," -").split()
          year=int(words[0])
          month=int(words[1])
          date=year+(month-0.5)/12
          lineno +=1
          if Year1:
            if month == 1:
              self.__firstFullYear = year
              Year1=False
              self.__lofull=lineno
          else:
            if month == 12:
              self.__lastFullYear=year
              self.__hifull=lineno
              
          self.__dates[date]=words[:2]
          self.__data.append(np.array(words[2:]))
        
    self.__data=np.array(self.__data,dtype=float)        
    self.__headers=Headers[2:]
    
  def headers(self):
    return self.__headers
    
  def monthly(self, column, month):
    if type(column)==type(int):
      col=column
    else:
      col=self.__headers.index(column)
    if month==None:
      return self.__data[self.__lofull-1:self.__hifull:,col]      
    else:
      return self.__data[self.__lofull+month-2:self.__hifull:12,col]

  def annual(self, column):
    if type(column)==type(int):
      col=column
    else:
      col=self.__headers.index(column)
    return np.mean([self.__data[self.__lofull+mo-1:self.__hifull:12,col] for mo in range(12)],axis=0)
    
  def data(self):
    return self.__data
      
  def filename(self):
    return self.__filename
    
if __name__ == "__main__":
  uah=UAH_MSU('tlt') 
  print uah.headers()
  print uah.data() 
  print uah.monthly('Globe',1)
  print uah.annual('Globe')

