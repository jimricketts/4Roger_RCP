# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 18:26:16 2015

@author: s4493222
"""


import os
import numpy as np
import regress
SVNRevision="$Revision: 371 $"
#doANCOVA
#http://r-eco-evo.blogspot.com.au/2011/08/comparing-two-regression-slopes-by.html


#set up the first anova to test model with break against one without, rquery saved as rqueryold
header= '''

data=read.csv("C:\\Users\\s4493222\\Documents\\abrupt\\SymmetricRegression_SVN_325\\Hadly2.csv")
outf <- file("C:\\Users\\s4493222\\Documents\\abrupt\\4Roger_TrendProbs\\result.ancova")
        '''
rqery=  '''
mod1<-aov(Anom~Year2*Class,data=data)
mod2<-aov(Anom~Year2+Class,data=data)
sink(outf)
data
summary(mod1)
summary(mod2)
anova(mod1,mod2)
priors <- subset(data, Year2 <=0)
posts <- subset(data, Year2 > 0)

reg1 <- lm(Anom~Year2,data=priors)
reg2 <- lm(Anom~Year2,data=posts)
summary(reg1)
summary(reg2)
reg.todo <- lm(Anom ~ Class*Year2 ,data=data)
summary(reg.todo)
mod3<-aov(Anom~Year2,data=data)
anova(mod1,mod3)
sink()
close(outf)
      '''
rqeryold=  '''
mod1<-aov(Anom~Year2*Class,data=data)
mod2<-aov(Anom~Year2+Class,data=data)
sink(outf)
data
summary(mod1)
summary(mod2)
anova(mod1,mod2)
priors <- subset(data, Year2 <=0)
posts <- subset(data, Year2 > 0)

reg1 <- lm(Anom~Year2,data=priors)
reg2 <- lm(Anom~Year2,data=posts)
summary(reg1)
summary(reg2)
reg.todo <- lm(Anom ~ Class*Year2 ,data=data)
summary(reg.todo)
sink()
close(outf)
      '''

sample= '''
            Df Sum Sq Mean Sq F value Pr(>F)    
Year2        1 2.0577  2.0577 266.809 <2e-16 ***
Class        1 0.0377  0.0377   4.884 0.0327 *  
Year2:Class  1 0.0247  0.0247   3.203 0.0809 .  
Residuals   41 0.3162  0.0077                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
            Df Sum Sq Mean Sq F value Pr(>F)    
Year2        1 2.0577  2.0577 253.511 <2e-16 ***
Class        1 0.0377  0.0377   4.641  0.037 *  
Residuals   42 0.3409  0.0081                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Analysis of Variance Table

Model 1: Anom ~ Year2 * Class
Model 2: Anom ~ Year2 + Class
  Res.Df     RSS Df Sum of Sq      F  Pr(>F)  
1     41 0.31621                              
2     42 0.34091 -1 -0.024703 3.2031 0.08089 .
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Call:
lm(formula = Anom ~ Year2, data = priors)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.19708 -0.07117 -0.01850  0.07484  0.14916 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept) 0.268732   0.038271   7.022 2.32e-07 ***
Year2       0.014745   0.002389   6.173 1.87e-06 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.09668 on 25 degrees of freedom
Multiple R-squared:  0.6038,	Adjusted R-squared:  0.588 
F-statistic:  38.1 on 1 and 25 DF,  p-value: 1.87e-06


Call:
lm(formula = Anom ~ Year2, data = posts)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.13361 -0.03542  0.01081  0.05131  0.12063 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept) 0.407754   0.032492  12.549 1.07e-09 ***
Year2       0.006617   0.003263   2.028   0.0595 .  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.07182 on 16 degrees of freedom
Multiple R-squared:  0.2045,	Adjusted R-squared:  0.1548 
F-statistic: 4.113 on 1 and 16 DF,  p-value: 0.05955


Call:
lm(formula = Anom ~ Year2 * Class, data = data)

Residuals:
      Min        1Q    Median        3Q       Max 
-0.197078 -0.069351  0.004777  0.065722  0.149158 

Coefficients:
                Estimate Std. Error t value Pr(>|t|)    
(Intercept)     0.407754   0.039731  10.263 6.81e-13 ***
Year2           0.006617   0.003990   1.659   0.1048    
ClassPre       -0.139022   0.052793  -2.633   0.0119 *  
Year2:ClassPre  0.008128   0.004542   1.790   0.0809 .  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.08782 on 41 degrees of freedom
Multiple R-squared:  0.8702,	Adjusted R-squared:  0.8607 
F-statistic: 91.63 on 3 and 41 DF,  p-value: < 2.2e-16
        '''
        
import statbreaks
import lowess
import numpy as np

tmpInputFile="ancova.in"
tmpQry="ancova.R"
tmpResults="ancova.out"

class ANCOVA(object):
  def __init__(self, Years, ys, breakyr,Rpath='"C:\\Program Files\\R\\R-3.0.2\\bin\\Rscript.exe"'):
  #def __init__(self, Years, ys, breakyr,Rpath='"C:\\Program Files\\R\\R-3.0.2\\bin\\x64\\Rscript.exe"'):
    self.__breakYear=breakyr
    self.__firstYear=min(Years)
    self.__lastYear=max(Years)
    #write a datafile
    xs1=[]
    xs2=[]
    ys1=[]
    ys2=[]
    with open(tmpInputFile, "w") as data:
      print >>data,"Year2,Anom,Class,Year"
      for i in range(len(Years)):
        y2=Years[i]-breakyr
        if y2 <=0.0:
          Class='Pre'
          xs1.append(y2)
          ys1.append(ys[i])
        else:
          Class='Post'
          xs2.append(y2)
          ys2.append(ys[i])
        print >>data,str(y2)+","+str(ys[i])+","+Class+","+str(Years[i])
    #compute regression stats my way first
   # self.__stats1=regress.analysed_regress(np.array(ys1), np.array(xs1))
   # self.__stats2=regress.analysed_regress(np.array(ys2), np.array(xs2))
    
    #compose the Rquery if we have to
    with open(tmpQry, "w") as qry:
      print >>qry,'data=read.csv("'+tmpInputFile+'")'
      print >>qry,'outf <- file("'+tmpResults+'")'
      print >>qry,rqery
    
    #execute the query
    query = Rpath + " "+tmpQry
    #print query
    try:
      os.system(query)
      with open(tmpResults, 'r') as res:
        self.__lines=res.readlines()
    except Exception as e:
      print str(e)
      
  def parseResult(self, lines):
    def tofloat(s):
      if s == "<":
        f=2e-16
      elif s == "NA":
        f=np.NaN
      else:
        try:
          f = float(s)
        except:
          if s[0] in ['<','>']:
            try:
              f = float(s[1:])
            except Exception as e:
              print "Cannot do tofloat(",s,")"
          pass
      return f
      
    state =-1
    nexttarget="2"
    for line in lines:
      if state == -1:
        if line[:3] == "Yea":
          state=0
      if state == 0:
        #significance of interaction with year
        if line[0]=='2':
          self.__prSignificantInteraction=float(line.split()[6])
          state = 1
          nexttarget= "(Intercept)"
      elif state == 1:
      #slopes and intercepts and their significance
        if line[:len(nexttarget)]==nexttarget:
          (_,self.__priorIntercept,_,_,self.__prPriorIntercept)=line.split()[:5]
          self.__priorIntercept=float(self.__priorIntercept)
          self.__prPriorIntercept=tofloat(self.__prPriorIntercept)
          state = 2
          nexttarget="Year2"
      elif state==2:
        if line[:len(nexttarget)]==nexttarget:
          (_,self.__priorSlope,_,_,self.__prPriorSlope)=line.split()[:5]
          self.__priorSlope=float(self.__priorSlope)
          self.__prPriorSlope=tofloat(self.__prPriorSlope)
          state = 2.5
          nexttarget="Multiple R-squared:"
      elif state == 2.5:
        if line[:len(nexttarget)]==nexttarget:
          self.__PriorRsq=float(line.split()[2][:-1])
          state = 3
          nexttarget="(Intercept)"
      elif state == 3:
        if line[:len(nexttarget)]==nexttarget:
          (_,self.__postIntercept,_,_,self.__prPostIntercept)=line.split()[:5]
          self.__postIntercept=float(self.__postIntercept)
          self.__prPostIntercept=tofloat(self.__prPostIntercept)
          state = 4
          nexttarget="Year2"
      elif state==4:
        if line[:len(nexttarget)]==nexttarget:
          (_,self.__postSlope,_,_,self.__prPostSlope)=line.split()[:5]
          self.__postSlope=float(self.__postSlope)
          self.__prPostSlope=tofloat(self.__prPostSlope)
          state = 4.5
          nexttarget="Multiple R-squared:"
      elif state == 4.5:
        if line[:len(nexttarget)]==nexttarget:
          self.__PostRsq=float(line.split()[2][:-1])
          state=5
          nexttarget="(Intercept)"
      #this will now be the ancova proper
      elif state == 5:
        if line[:len(nexttarget)]==nexttarget:
          (_,self.__totalShift,_,_,self.__prTotalShift)=line.split()[:5]
          self.__totalShift=float(self.__totalShift)
          self.__prTotalShift=tofloat(self.__prTotalShift)
          state = 6
          nexttarget="ClassPre"
      elif state == 6:
        if line[:len(nexttarget)]==nexttarget:
          (_,self.__diffShift,_,_,self.__prShift)=line.split()[:5]
          self.__diffShift=-float(self.__diffShift)
          self.__prDiffShift=tofloat(self.__prShift) #this is Pr intercepts differ
          self.__priorIntercept=self.__totalShift-self.__diffShift
          state = 7
          nexttarget="Year2"
      elif state == 7:
        if line[:len(nexttarget)]==nexttarget:
          (_,self.__postTrend,_,_,self.__prPostTrend)=line.split()[:5]
          self.__postTrend=float(self.__postTrend)
          self.__prPostTrend=tofloat(self.__prPostTrend)
          state = 8
          nexttarget="ClassPre:Year2"
      elif state == 8:
        if line[:len(nexttarget)]==nexttarget:
          (_,self.__diffTrend,_,_,self.__prDiffTrend)=line.split()[:5]
          self.__diffTrend=float(self.__diffTrend)
          self.__prDiffTrend=tofloat(self.__prDiffTrend)
          self.__PreTrend=self.__prPostTrend + self.__prDiffTrend
          nexttarget="2"
          state = 9
      elif state == 9: #using an anova between break model and linear model
        if line[:len(nexttarget)]==nexttarget:
          self.__prbreak=line.split()[6]
          self.__prbreak=float(self.__prbreak)
          state = 10
          
    return state         
    
  def formatted(self):
    lines=["Summary of a break year after "+str(self.__breakYear)]
    lines.append("Priors")
    lines.append("Slope :%01.4f (Pr %g), Intercept at Break :%01.4f (Pr %g), Mult-Rqs:%g" % (self.__priorSlope,self.__prPriorSlope,self.__priorIntercept,self.__prPriorIntercept,self.__PriorRsq))
#    lines.append("Slope :%01.4f (Pr %g), Intercept at Break :%01.4f (Pr %g), Mult-Rqs:%g OLS-Rsq: %g, Pr(%g)" % (self.__priorSlope,self.__prPriorSlope,self.__priorIntercept,self.__prPriorIntercept,self.__PriorRsq,self.__stats1['rsq'],self.__stats1['prob']))
    lines.append("Posts")
    lines.append("Slope :%01.4f (Pr %g), Intercept at Break :%01.4f (Pr %g), Mult-Rsq:%g" % (self.__postSlope,self.__prPostSlope,self.__postIntercept,self.__prPostIntercept,self.__PostRsq))
#    lines.append("Slope :%01.4f (Pr %g), Intercept at Break :%01.4f (Pr %g), Mult-Rsq:%g OLS-Rsq: %g, Pr(%g)" % (self.__postSlope,self.__prPostSlope,self.__postIntercept,self.__prPostIntercept,self.__PostRsq, self.__stats2['rsq'],self.__stats2['prob']))
    lines.append("Change Point")
    lines.append("Change of Slope: %g Pr(%g) " % (-self.__diffTrend, self.__prDiffTrend))
    lines.append("SHIFT at  break: %g Pr(%g) " % (self.__diffShift, self.__prDiffShift))
    lines.append("Pr of No Break : %g " % (self.__prbreak, ))

    return lines
    
  def results(self):
    return self.__lines
    
  def firstYear(self):
    return self.__firstYear
    
  def lastYear(self):
    return self.__lastYear

  def firstSlope(self):
    try:
      return self.__priorSlope
    except:
      print "o0ps"
    
  def firstSlopeProb(self):
    return self.__prPriorSlope

  def secondSlope(self):
    return self.__postSlope
    
  def secondSlopeProb(self):
    return self.__prPostSlope
    
  def changeOfTrend(self):
    return self.__diffTrend
    
  def changeOfTrendProb(self):
    return self.__prDiffTrend
    
  def shift(self):
    return self.__diffShift
    
  def trendProb(self):
    return self.__prDiffTrend
    
  def firstRsq(self):
    return self.__PriorRsq

  def secondRsq(self):
    return self.__PostRsq
    
  def prBreak(self):
    return self__prbreak

if __name__ == "__main__":

  fn1=os.environ["HOMEPATH"]+"\\Documents\\abrupt\\c_test\\GISSTEMP_v3_April_2015_withDiffs_from_2014_24N44N\\GISSTEMP_v3_April_2015_withDiffs_from_2014.csv_24N44N_0.trace"
  #os.environ["HOMEPATH"]+"\\Documents\\abrupt\\4Roger_Nature_SVN_264\\HadCRUT.4.2.0.0.annual_ns_avg//HadCRUT.4.2.0.0.annual_ns_avg.txt_0.trace")    
  bp=statbreaks.brkrpt(fn1)
  yearsegs, ysegs, xsegs = bp.segments()
  pre98Years=yearsegs[-2]
  pre98temps=ysegs[-2]
  post98Years=yearsegs[-1]
  post98temps=ysegs[-1]
  Years=list(pre98Years[:])
  Years.extend(list(post98Years[:]))
  ys=list(pre98temps[:])
  ys.extend(list(post98temps[:]))
  
  print bp.breaks()
  
  ancova=ANCOVA(Years, ys, bp.breaks()[-2], Rpath='"C:\\Program Files\\R\\R-3.0.2\\bin\\Rscript.exe"')
  print ancova.parseResult(ancova.results())
  for l in ancova.formatted():
    print l
  