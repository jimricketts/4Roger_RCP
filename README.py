# -*- coding: utf-8 -*-
"""
Created on Mon Mar 09 12:49:39 2015

@author: s4493222
"""
SVNRevision="$Revision: 308 $"
text= '''
This is a readme file for the b_test directory.
Two issues.
1) Does the recursive test fail on the right too readily? YES!!
2) Is the instability in GISSTEMPtp62-13v_64S44S real or an issue?  NO!! An ISSUE - fixed

Experiment 1. The issue here is that a number of 21st C runs show a "final" break in the recusive test that the 
convergent test finds after.

A problem exists

  #this is aexperimet 1. Load a
  fn="I:\\4Roger_ClimatePapers_SVN_280_tas_Amon\\tas_Amon_CESM1-CAM5_rcp60_r1i1p1_200601-210012\\tas_Amon_CESM1-CAM5_rcp60_r1i1p1_200601-210012.GW_0.trace"
  sb=statbreaks.brkrpt(fn)
  years=sb.years()
  ys=sb.ys()
  xs=sb.xs()
  breaks=sb.breaks()
  recbr=sb.recursive()
  print recbr, breaks
  
  bv=bivariate.bivariate(ys, xs, anomalise=False)
  bv1=bivariate.bivariate(ys[bv.maxIndexTi()+1:], xs[bv.maxIndexTi()+1:], anomalise=False)
  [1850.0, 1883.0, 1914.0, 1969.0, 1997.0, 2026.0, 2100.0] [1850.0, 1883, 1914, 1971, 1997, 2013, 2025, 2037, 2053, 2062, 2075, 2085, 2100.0]
>>> bv1.maxTi()
56.181595670461455
>>> years[bv.maxIndexTi()+1:]
array([ 2027.,  2028.,  2029.,  2030.,  2031.,  2032.,  2033.,  2034.,
        2035.,  2036.,  2037.,  2038.,  2039.,  2040.,  2041.,  2042.,
        2043.,  2044.,  2045.,  2046.,  2047.,  2048.,  2049.,  2050.,
        2051.,  2052.,  2053.,  2054.,  2055.,  2056.,  2057.,  2058.,
        2059.,  2060.,  2061.,  2062.,  2063.,  2064.,  2065.,  2066.,
        2067.,  2068.,  2069.,  2070.,  2071.,  2072.,  2073.,  2074.,
        2075.,  2076.,  2077.,  2078.,  2079.,  2080.,  2081.,  2082.,
        2083.,  2084.,  2085.,  2086.,  2087.,  2088.,  2089.,  2090.,
        2091.,  2092.,  2093.,  2094.,  2095.,  2096.,  2097.,  2098.,
        2099.,  2100.])
>>> years[bv.maxIndexTi()+1:][bv1.maxIndexTi()]
2069.0   <<This should be in the recursive set
      
      '''
  
def getreturn(astring):
  pos1=astring.find("-> ")
  pos2=astring.find("] [")
  print astring[pos1+3:pos2+1]
    
if __name__ == "__main__":
  import ConvergentBreaks as convergent_breaks
  import bivariate_multi as bivariate
  import statbreaks
  import random
  
  with open ("giss.x", "r") as f:
    lines = f.readlines()
    for line in lines:
      getreturn(line)
  
  
  #print text
  '''
  #this is aexperimet 1. Load a
  fn="I:\\4Roger_ClimatePapers_SVN_280_tas_Amon\\tas_Amon_CESM1-CAM5_rcp60_r1i1p1_200601-210012\\tas_Amon_CESM1-CAM5_rcp60_r1i1p1_200601-210012.GW_0.trace"
  fn="C:\\Users\\s4493222\\Documents\\abrupt\\4Roger_ClimatePapers_SVN_280\\GISSTEMPto6-2013b_64S44S\\GISSTEMPto6-2013b.csv_64S44S_0.trace"  
  sb=statbreaks.brkrpt(fn)
  years=sb.years()
  ys=sb.ys()
  xs=sb.xs()
  breaks=sb.breaks()
  recbr=sb.recursive()
  print recbr, breaks
  
  bv=bivariate.bivariate(ys, xs, anomalise=False)
  print bv.allPoints(pr=0.01)
  #bv1=bivariate.bivariate(ys[bv.maxIndexTi()+1:], xs[bv.maxIndexTi()+1:], anomalise=False)
  for i in range(100):
    print "Iteration", i
    xs = np.array([random.random() for y in ys])
    convergent_breaks.TraceFile = "test_b_GISSTEMPto6-2013b_64S44S_"+str(i)+'.trace'
    cb=convergent_breaks.convergentBreaks(ys, xs, years, 'GWAnom', mode="control", guide="Stability",screenpr=0.01, pr=0.01, trace=True)
  
  
  "tas_Amon_CESM1-CAM5_rcp60_r1i1p1_200601-210012.GW_99.trace:Returning [1850.0, 1914.0, 1969.0, 1997.0, 2026.0, 2100.0] -> [1850.0, 1883, 1914, 1971, 1997, 2013, 2025, 2037, 2053, 2068, 2100.0] [((1883.0, 0.0), (14.875033, 0.64558715), (-0.11563264210073956, 0.0027004624096281888)), ((1914.0, 0.0), (49.033031, 0.78020418), (0.21243845599012101, 0.0027367499834629015)), ((1971.0, 0.0), (53.245468, 0.54289806), (0.23785867534437785, 0.0019239136220113458)), ((1997.0, 0.0), (33.447506, 0.28639904), (0.44813400369110218, 0.006210219345962842)), ((2013.0, 0.0), (19.401585, 0.50304854), (0.35545022591195463, 0.010183146070949856)), ((2025.0, 0.0), (15.055582, 0.54039711), (0.25107228004410143, 0.0061754873746148811)), ((2037.0, 0.0), (17.399935, 0.3711513), (0.32611510797413823, 0.0075495735110717808)), ((2053.0, 0.0), (16.207022, 0.48892793), (0.31542020438022467, 0.0085391515289537819)), ((2068.0, 0.0), (17.899317, 0.56454796), (0.51371517074845496, 0.016539671961291948))]"

  f="C:\\Users\\s4493222\\Documents\\abrupt\\4Roger_ClimatePapers_SVN_280\\GISSTEMPto6-2013b_64S44S\\GISSTEMPto6-2013b.csv_64S44S_0.trace"
  ''' 
