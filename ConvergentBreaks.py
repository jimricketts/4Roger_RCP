# -*- coding: utf-8 -*-
"""
Created on Wed Oct 01 09:32:41 2014
ConvergentBreaks take 2
This time with a proper merge
@author: s4493222
"""


import numpy as np
import bivariate_multi as bivariate
import scipy.stats as stats
import random
import recursetest
import copy
import regress

from bisect import bisect_left
import AIC
SVNRevision="$Revision: 317 $"

TraceFile = ""
#import shuffle#
def takeClosest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    from : http://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return after
    else:
       return before

def takeClosestIndex(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    from : http://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return pos
    if pos == len(myList):
        return pos-1
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return pos
    else:
       return pos-1
       
class CBException(Exception):
  pass
  
  def __str__(self):
    return "CBException: "+str(self.msg)


def resample_break(testdata, datayears, N=30, withmode=False):
  brks=[]
  tis=[]
  shifts=[]
  try:
    for i in range(N):
      step=0
      controldata = np.array([random.random() for y in datayears])
      step=1
      bv=bivariate.bivariate(testdata,controldata, anomalise=False, pr=0.01)
      step=2
      tis.append(bv.maxTi())
      step=3
      brks.append(datayears[bv.maxIndexTi()])
      step=4
      shifts.append(bv.stepChange())
      step=5
  except Exception as e:
    print str(e), step
    raise e
  if withmode:
    yearfreqs=np.bincount(brks)
    #print brks,len(yearfreqs),yearfreqs,range(int(datayears[0])-1,int(datayears[-1])-1)
    first=second=firstval=secondval=0
    for i in range(int(min(brks))-1,int(max(brks))+1):
     # print i
      if yearfreqs[i] >firstval:
        second=first
        first=i
        firstval=yearfreqs[first]
        secondval=yearfreqs[second]
      elif yearfreqs[i] >secondval:
        second=i
        secondval=yearfreqs[second]
    #now we also will need the mean of the Tis and shifts for each of the modal values
    brks=np.array(brks)
    mask1=np.where(brks==first) 
    timean1=np.mean(np.array(tis)[mask1])
    shmean1=np.mean(np.array(shifts)[mask1])
    if second == 0:
      return stats.norm.fit(brks[mask1]), stats.norm.fit(np.array(tis)[mask1]),stats.norm.fit(np.array(shifts)[mask1]), [(first, float(firstval)/N), None, None,(timean1, shmean1), None]
    else:
      mask2=np.where(brks==second) 
      timean2=np.mean(np.array(tis)[mask2])
      shmean2=np.mean(np.array(shifts)[mask2])
      return stats.norm.fit(brks[mask1]), stats.norm.fit(np.array(tis)[mask1]),stats.norm.fit(np.array(shifts)[mask1]), [(first, float(firstval)/float(N)), (second, float(secondval)/float(N)), yearfreqs[-len(datayears):],(timean1, shmean1), (timean2, shmean2)]
#    if second == 0:
#      return stats.norm.fit(brks), stats.norm.fit(tis),stats.norm.fit(shifts), [(first, float(firstval)/N), None, None,(timean1, shmean1), None]
#    else:
#      mask2=np.where(brks==second) 
#      timean2=np.mean(np.array(tis)[mask2])
#      shmean2=np.mean(np.array(shifts)[mask2])
#      return stats.norm.fit(brks), stats.norm.fit(tis),stats.norm.fit(shifts), [(first, float(firstval)/float(N)), (second, float(secondval)/float(N)), yearfreqs[-len(datayears):],(timean1, shmean1), (timean2, shmean2)]
  else:
    return stats.norm.fit(brks), stats.norm.fit(tis),stats.norm.fit(shifts)  
  
  
def convergentBreaks_Inner(testdata, controldata, datayears, aicControl, model, trace=True, shallow=True, keepFirst=False, pr=0.01, screenpr=0.05, terminalsolution=True, minInterval=7):
  '''
  This code attempst to iteratively test all 
  '''
  #The issue with convergence is that whether a breakpoint is admitted to the yearly breaks is determined 
  #initial analysis over full data
  #implement as (yet another) statemachine
  
  def statelabel(state1, state2):
    return "state["+str(state1)+"]->["+str(state2)+"]"
  
  
  
  def merge_pass(currentbreaks, testdata, controldata, datayears):
    prevbreaks=copy.copy(currentbreaks)
    newlist=[prevbreaks.pop(0)]    
    lo=takeClosestIndex(datayears, newlist[0])
    print >>tf, "considering",str(prevbreaks)
    #JHR 13/2/15 Here we consider a prolog and an epilog since it ius necessary to test within the 
    #first and last spans, just to see if something is shielded.
    PrologDone=False
    
    try:
      if len(currentbreaks) < 2:
        print "merge_pass.currentbreaks=",currentbreaks
      if len(prevbreaks) > 1:
        hi = takeClosestIndex(datayears, prevbreaks[1])+1
      else: #JHR SVN 262 15/1/2015 If called with only the end data s bounds this will crash
        hi = takeClosestIndex(datayears, prevbreaks[0])+1
        PrologDone = True #so just prevent doing the same interval twice should we examine a provionally empty span
    except:
      print "merge_pass.prevbreaks=",prevbreaks,"merge_pass.currentbreaks=",currentbreaks
      raise
    #SVN 280 JHR 13/2/15 Consider
    firsthi = takeClosestIndex(datayears, prevbreaks[0])+1

    statlist = []
    state=0
    safetystep = 0
    Terminating = False #Flag to say we are on last pass
    while state >= 0: #state < 0 indicates termination
      
      #state evaluation loop is always get new list, then depending on state process
      if state==0: #ythen we have commenced
        #It should not happen, but was getting exceptions causing hi to be wrongin termination case
        #JHR 2/3/2015
        if Terminating:
          hi = takeClosestIndex(datayears, currentbreaks[-1])+1
        rt1=recursetest.recurse(testdata[lo:hi], controldata[lo:hi], datayears[lo:hi], model, smooth=False, trim = 0, pr=screenpr, anom=False, withshifts=True)
        candidates= np.sort(rt1.breakyears().keys()).tolist() 
        if Terminating and candidates != [] and trace:
          print >>tf, "Final pass found ",candidates, " between ", datayears[[lo,hi-1]], "state[0] -> state[1]"  
        if not PrologDone:
          PrologDone = True #so revert to previous code and commence testing spans
          rtprolog=recursetest.recurse(testdata[lo:firsthi], controldata[lo:firsthi], datayears[lo:firsthi], model, smooth=False, trim = 0, pr=screenpr, anom=False, withshifts=True)
          firstcandidates= np.sort(rtprolog.breakyears().keys()).tolist() 
          if firstcandidates != []:
            if trace:
              print >>tf, "B4 found ",firstcandidates, " between ", datayears[[lo,firsthi-1]], "state[0] -> state[0]"  
            prevbreaks.insert(0,firstcandidates[0])
            hi = firsthi
            state = 0
            #The following can cause deletion of a year, better to just start over 
            #firstcandidates.extend(candidates)
            #candidates=np.sort(list(set(firstcandidates))).tolist()
          else:
            state = 1
        else:
          state = 1 

      elif state == 1: #the state after we have a new list, candidates needs to be evaluated, firstly how many?
        if len(candidates) == 0:
          state = 2 #the state of do a drop 
        elif len(candidates) == 1:
          state = 3 # the state of process just one
        else:
          state = 4 # the state of process multiples

      elif state == 2: #then our test segment yielded no acceptable break and we will just drop it
        if trace:
          if len(prevbreaks) > 1:
            print >>tf, "** Break ",prevbreaks[1],"no longer between ",datayears[lo], datayears[hi-1], statelabel(state,5)
          else:
            print >>tf, "** Break ",prevbreaks,"no longer between ",datayears[lo], datayears[hi-1], statelabel(state,5)
        state = 5 #the state of a break was dropped by itself (therefore we don't move lo)

      elif state == 3:# the state of process just one
        #so now we set up an evaluation of one element.
        #we rely on candidates, lo, and hi
        modes=[None,None]
        ystats, tstats, shiftstats,modes = resample_break(testdata[lo:hi], datayears[lo:hi],N=100,withmode=True)
        #ystats, tstats, shiftstats = resample_break(testdata[lo:hi], datayears[lo:hi],withmode=False)
        #now the tstats may mean it is no longer significant
        #or the ystats may indicate it's too close to the previous
        #or there may be another issue unthought about
        #how about the Ti0?
        
        crit30=bivariate.critTi(pr, max(minInterval,min(30, 1+datayears[hi-1]-datayears[lo])))
        #now the retested point may have moved
        if modes[-2][0] < crit30: #then it is not significant (was tstats[0])
          if trace:
            print >>tf,  "-- Candidate ", candidates[0], "(",modes[0][0],")(",modes[-2][0],") did not exceed ", crit30, "between ",datayears[lo], datayears[hi-1], "state[",state,"]" ,
          state= 6 #the state of a break was dropped because it was a problem
          if trace:
            print >>tf, "-> state[",state,"]"        
        elif modes[0][0] -  (datayears[lo] -1) < minInterval:#JHR This was one year out, last break is datayears[0] -1 #was ystat[0] then it probably moved and in any case is too close to last break
          if trace:
            print >>tf,  "<< Candidate ", candidates[0], "(",modes[0][0],")(",modes[-2][0],") too close to ", datayears[lo], "between ",datayears[lo], datayears[hi-1] , "state[",state,"]",
          state= 5.5 #the state of two peaks close togeter - choose 1
          if trace:
            print >>tf, "-> state[",state,"]"        
        # hold this thought
        #elif (ystats[0] - 2 * ystats[1] < min(testyr, datayears[lo]) or ystats[0] + 2 * ystats[1] > max(testyr, datayears[hi-1])):
        #now we consider a bunch of rules about the modality if we have done this
        elif modes[1] != None: #so there's more than one choice
          #if the first mode is more than 90% (Roger's verbal rule) the all is OK
          if modes[0][1] >= 0.9:
            state = 7 #the state of save the stats and then move on
          elif modes[0][1]+modes[1][1] < 0.7: #if less than 70% of the values are between two modes then it must be blurry
            if trace:
              print >>tf,  "~~ Candidate no strong modes ", candidates[0], "(",modes,") less than 70% ", datayears[lo], "between ",datayears[lo], datayears[hi-1] , "state[",state,"]",
            state = 6  #the state of a break was dropped because it was a problem
            if trace:
              print >>tf, "-> state[",state,"]"        
          elif modes[1][1] <= 0.2: #so at half the time the mode is one value
            state = 7 #the state of save the stats and then move on
          elif abs(modes[0][0] - modes[1][0]) <= minInterval /2.0:
            if trace:
              print >>tf,  "~+ Candidate accepted with two modes ", candidates[0], "(",modes,") < minInterval /2.0 ", datayears[lo], "between ",datayears[lo], datayears[hi-1] , "state[",state,"]",
            state = 7
            if trace:
              print >>tf, "-> state[",state,"]"        
          elif abs(modes[0][0] - modes[1][0]) > minInterval:
            #so for now we will override the candidate list and treat as two candidates _NOT ANY MORE
            if trace:
              print >>tf,  "~~ Candidate two strong modes - treat as both possible ", candidates[0], "(",modes,") > ", datayears[lo], "between ",datayears[lo], datayears[hi-1] , "state[",state,"]",
            #candidates = [min(modes[0][0], modes[1][0])]#, max(modes[0][0], modes[1][0])]
            state = 9.5                  
            if trace:
              print >>tf, "-> state[",state,"]"        
          else:
            if trace:
              print >>tf,  "~~ Candidate two strong but close modes ", candidates[0], "(",modes,") > minInterval /2.0 ", datayears[lo], "between ",datayears[lo], datayears[hi-1] , "state[",state,"]",
            state = 7 #This was a definite error - had been dropping these. setting to state 6
            if trace:
              print >>tf, "-> state[",state,"]"        
          #that's basically it for mode processing
          #however one more thing to check              
        if state == 3: #final pass through
          if abs(candidates[0] - ystats[0]) > 2.0 * ystats[1]:#resampling shows bad point JHR: 20150107 BAD TEST FIXME!
            if trace:
              print >>tf,  "<<! Candidate ", candidates[0], "(",ystats,")(",tstats,") too far from ", candidates[0], "between ",datayears[lo], datayears[hi-1] , ". Will replace, state[",state,"]",
            state= 9 #the state of a break was changed due to resampleing
            if trace:
              print >>tf, "-> state[",state,"]"        
          else:
            state = 7 #the state of save the stats and then move on
        if state == 7 and trace:
            if candidates[0] != ystats[0]:
              print >>tf, "Candidate has moved : modes are " ,modes
            print >>tf,  "=3! Candidate ", candidates[0], "(",ystats,")(",tstats,") OK between ",datayears[lo], datayears[hi-1] , " state[3] - state[7]"
      elif state == 4:
        #test the second candidate to see is it's going to be stable. 
        #But we first need to know what the first candidate is.
      
      
        modes1=[None,None]
        lo1 = takeClosestIndex(datayears, candidates[0])+1        
        hi1 = takeClosestIndex(datayears, candidates[1])+1
        #SO WHAT IS the new lo bound?
        #If this one is above critical use it, otherwise don't. Use the second candidate which by virtue of how is was generated is likely to persist
        ystats0, tstats0, shiftstats0,modes0 = resample_break(testdata[lo:hi1], datayears[lo:hi1],N=100,withmode=True)
        lo0 = takeClosestIndex(datayears, modes0[0][0])+1
        #so now we use the new low bound against the proposed upper bound 
        ystats1, tstats1, shiftstats1,modes1 = resample_break(testdata[lo0:hi], datayears[lo0:hi],N=100,withmode=True)
        crit30a=bivariate.critTi(pr, max(minInterval,min(30, 1+datayears[hi-1]-datayears[lo0])))
        lowcandidate = candidates[0]
        if modes1[-2][0] >= crit30a: #then we save that one to use
          if trace:
            print >>tf,  "!-- Candidate ", candidates[0], " reset to ",modes0[0][0], " after resample gave ",modes0
          lowcandidate = modes0[0][0]
          lo1 = lo0
        else: #otherwise need to make a decision with the unchnaged candidate 
          ystats1, tstats1, shiftstats1,modes1 = resample_break(testdata[lo1:hi], datayears[lo1:hi],N=100,withmode=True)
          crit30a=bivariate.critTi(pr, max(minInterval,min(30, 1+datayears[hi-1]-datayears[lo1])))

        if modes1[-2][0] < crit30a: #then it is not significant (was tstats[0])
          if trace:
            print >>tf,  "!-- Candidate ", candidates[1], "(",modes1[0][0],")(",modes1[-2][0],") did not exceed ", crit30a, "between ",datayears[lo1], datayears[hi-1], " will not use: state[",state,"]" ,
        else:
          #insert       
          if trace and candidates[1] in prevbreaks:
            print >>tf, "KNOWN DUPLICATE GOING INTO LIST", candidates[1], prevbreaks, "THIS IS OK"  
          if trace:
            print >>tf,  "++ (new) Candidate ", candidates[0], "will be tested and ", candidates[1],"inserted", prevbreaks, "becomes", 
          prevbreaks.insert(1, candidates[1]) #avoid relooping
          if trace:
            print >>tf, prevbreaks, "between ",datayears[lo], datayears[hi-1], "state[",state,"]",
          #lo=takeClosestIndex(datayears, newlist[-1])+1
          if not Terminating: hi = takeClosestIndex(datayears, prevbreaks[1])+1
          
        candidates=[lowcandidate]
        state = 3
        if trace:
          print >>tf, "-> state[",state,"]"        

#      elif state == 4.1: #Old code
#        if trace:
#          print >>tf,  "++ (new) Candidate ", candidates[0], "will be tested and ", candidates[1],"inserted", prevbreaks, "becomes", 
#        prevbreaks.insert(1, candidates[1]) #avoid relooping
#        if trace:
#          print >>tf, prevbreaks, "between ",datayears[lo], datayears[hi-1], "state[",state,"]",
#        #lo=takeClosestIndex(datayears, newlist[-1])+1
#        hi = takeClosestIndex(datayears, prevbreaks[1])+1
#        candidates=[candidates[0]]
#        state = 3
#        if trace:
#          print >>tf, "-> state[",state,"]"        
#
#      elif state == 4.2: #working on
#        if trace:
#          print >>tf,  "++ (new) Candidate ", candidates[0], "will be tested and ", candidates[1],"inserted", prevbreaks, "becomes", 
#        prevbreaks[0]= candidates[0]
#        prevbreaks[1]= candidates[1] #avoid relooping
#        if trace:
#          print >>tf, prevbreaks, "between ",datayears[lo], datayears[hi-1], "state[",state,"]",
#        #lo=takeClosestIndex(datayears, newlist[-1])+1
#        hi = takeClosestIndex(datayears, prevbreaks[1])+1
#        state = 0
#        if trace:
#          print >>tf, "-> state[",state,"]"        
#

      elif state == 5:
        #the state of a break was dropped by itself (therefore we don't move lo)
        #actually we have to move it a bit to stop a loop if we then get it back - 
        #lo=takeClosestIndex(datayears, newlist[0])
        #do the dropping 
        if Terminating: 
          state = 10 # so avoid this processing
        else:
          prb=prevbreaks.pop(0)
          if len(prevbreaks) >1:
            hi = takeClosestIndex(datayears, prevbreaks[1])+1
            lo += minInterval/2
            if trace:
              old = hi-1
              print >>tf,  "!! Hi bound, ",datayears[old], " to give ",datayears[lo], datayears[hi-1], "state[",state,"]",
            state = 0
          else:
            if trace:
              print >>tf,  "!! Hi bound, ",datayears[hi-1], " last between",datayears[lo], datayears[hi-1], "state[",state,"]",
            state = 8 #nothing to do after dropping an upper bound
          if trace:
            print >>tf, "-> state[",state,"]"        

      elif state == 5.5:
        #the state of two peaks close togeter - choose 1 - if we choose the older one just drop this one
        #if we choose the later 1 delete reference to the prior one
        #similar code to state 3
        #but to prevent a loop we must prevent setting up the same interval again
        
        if len(statlist) == 0:
          state = 7 #accept it
          if trace: print >>tf,  "--!! Candidate ", candidates[0], "must be too close to start, but accept it, state[5.5] -> state[7]"
        else:
          if trace: print >>tf, "In state[5.5]"
          lo1=takeClosestIndex(datayears, newlist[-2])+1 #go back one step
          ystats, tstats, shiftstats,nmodes = resample_break(testdata[lo1:hi], datayears[lo1:hi], withmode=True)
          
          crit30=bivariate.critTi(pr, max(minInterval,min(30, 1+datayears[hi-1]-datayears[lo1])))
          #now the retested point may have moved
          if tstats[0] < crit30: #then it is not significant and this sia serious issue since one of these points previously was
            if trace:
              print >>tf,  "--!! Candidate ", candidates[0], "(",ystats,")(",tstats,") did not exceed ", crit30, "between ",datayears[lo], datayears[hi-1], "state[",state,"]" ,          
            state= 5 #the state of a break was dropped because it was a problem
            if trace:
              print >>tf, "->!! state[",state,"]"        
          else:
            if trace:
              
              print >>tf, "--!! go back to old break ",newlist[-2], "from ", datayears[lo], " candidate ",candidates[0], "new ", ystats[0], "between ",datayears[lo1], datayears[hi-1],"state[5.5]->",
            candidates=[round(ystats[0])]
            if ystats[0] -  datayears[lo] < minInterval:
              #we haven't moved very far, so set up a safety step
              safetystep = int(abs(ystats[0] -  datayears[lo]))
            if len(statlist) > 0: #if we have already got  
              try:
                newlist.pop(-1)
                statlist.pop(-1)
              #remove old dates
                state = 7 #the state of save the stats and then move on
              except: #we must not be able to pop, just delete it
                print >>tf, "no pop state [5.5] -> "
                state = 5
                pass
            else:
              state = 8 #nothing to do after dropping a candiate
            if trace: print >>tf, "state[",state,"]"
        #do the dropping 
        #alternative strategy when the new point is too close to the previous, go back and choose the 
        # highest Ti) of the two
      
     
      elif state == 6: #the state of a break was dropped because it was a problem so currently just update the indices past the identified candidate break
        #relies on candidate[0]
        if Terminating:
          state = 10
        else:
          oldlo = lo
          lo = takeClosestIndex(datayears, candidates[0])+1
          if len(prevbreaks) <2: #then cannot go further, we are done. JHR 15/1/2014 SVN262   
            if trace:
              print >>tf,  "$$ Lo bound ONLY updated to give ",datayears[lo], datayears[hi-1], "state[",state,"]",
            state = 0 #nothing to do after dropping an upper bound
            
          else:
            #need to see if we can move the upper bound safely retaining the low
            if len(prevbreaks)>2:
              hi3 = takeClosestIndex(datayears, prevbreaks[2])+1
            else:
              hi3 = takeClosestIndex(datayears, prevbreaks[1])+1
            ystats3, tstats3, shiftstats3,modes3 = resample_break(testdata[oldlo:hi3], datayears[oldlo:hi3],N=100,withmode=True)
            if modes3[1] == None or modes3[0][1]+modes3[1][1] >= 0.7:
              hi = hi3
              lo = oldlo
              state = 0 #no more to do but collect another point
              prb=prevbreaks.pop(0)
              if len(prevbreaks) >1 and prb == prevbreaks[0]: #JHR 17mar2015 was looping here.
                 prb=prevbreaks.pop(0) #remove the second instance if present to prevent a loop

              if trace:
                print >>tf,  "$$ Moved high bound, NOW ",datayears[oldlo] , datayears[hi-1], "dropped ", prb, "to give", prevbreaks, "state[",state,"]",
            
            #if this now violates the length rule update hi as well
            elif datayears[hi-1] - datayears[lo] < minInterval:
              if trace:
                print >>tf,  "$$ Lo bound, ",datayears[oldlo],"too close to Hi to use, cannot apply ",datayears[lo], datayears[hi-1], "state[",state,"]",
              state = 5#the state of a break was dropped by itself (therefore we don't move lo)
            else:
              if trace:
                print >>tf,  "$$ Lo bound, updated to give ",datayears[lo], datayears[hi-1], "state[",state,"]",
              state = 8 #nothing to do after dropping an upper bound
          if trace:
            print >>tf, "-> state[",state,"]"        

      elif state == 7:#the state of save the stats and then move on
        #we rely on candidates, lo, and hi, ystats, tstats, shiftstats
        if ystats == None:
          print "Trouble - this is a debug line "
        if not (Terminating and abs(currentbreaks[-1] -ystats[0]) <  minInterval):  
          statlist.append((ystats, tstats, shiftstats))
          newlist.append(int(round(ystats[0])))
        #newlist.append(candidates[0]) #Save the actual candidate
        if not Terminating: 
          lo = takeClosestIndex(datayears, newlist[-1])+1+safetystep
          if lo > len(datayears): lo -= safetystep          
          safetystep =0
          prb=prevbreaks.pop(0)
          if len(prevbreaks) >1 and prb == prevbreaks[0]:
            if trace:
              print >>tf, "Investigate : just popped", prb, prevbreaks, "state[7]"
            prb=prevbreaks.pop(0)  #JHR remove the duplicate as well. 10/2/2015
            
          if len(prevbreaks) >1:
            hi = takeClosestIndex(datayears, prevbreaks[1])+1
            state = 0 #round again
          else:
            state = 10 #nothing more to do
        else:
          state = 10
        if trace:
          print >>tf,  "\\\\", "Saved ",newlist[-1],statlist[-1], "state[7] -> state[",state,"]"
      
      elif state == 8: #nothing to do after dropping an upper bound
        if trace:
          print >>tf, "-- Update upper bound", datayears[hi-1], "to ",
        try:
          prb=prevbreaks.pop(0)
        except:
          if trace:
            print >>tf, "empty list",
          pass
        if not Terminating and len(prevbreaks) >1:
          hi = takeClosestIndex(datayears, prevbreaks[1])+1
          if trace:
            print >>tf,datayears[hi-1], "state[",state,"]",        
          state = 0 #round again
        else:
          state = 10 #nothing more to do
        if trace:
          print >>tf, "->","state[",state,"]"        
          
      elif state== 9: #the state of a break was changed due to resampleing
        if trace:
          print >>tf,  "&& candidate, ",candidates[0], " no longer updated to give ",
        #candidates[0] = round(ystats[0])
        if trace:
          print >>tf,round(ystats[0]),"state[",state,"]",
        state = 7
        if trace:
          print >>tf, "-> state[",state,"]"        

      elif state== 9.5: #the state of two strong modes record and exit
        if trace:
          print >>tf,  "&&-> candidate, ",candidates[0], "updated to give ",
        if modes[3][0] >= crit30:
          ystats=(modes[0][0],0.0)
          tstats=(modes[3][0], 0.0)
          shiftstats=(modes[3][1],0.0)
          candidates = [modes[0][0]]
          if trace:
            print >>tf, " FIRST MODE ",
        else:
          ystats=(modes[1][0],0.0)
          tstats=(modes[4][0], 0.0)
          shiftstats=(modes[4][1],0.0)
          candidates = [modes[1][0]]
          if trace:
            print >>tf, " SECOND MODE ",
        
        if trace:
          print >>tf,candidates[0],"state[",state,"]",
        state = 7
        if trace:
          print >>tf, "-> state[",state,"]"        

      elif state == 10:
        if Terminating:
          newlist.append(currentbreaks[-1])
          state = -1
        else:
          Terminating = True
          #lo = takeClosestIndex(datayears, newlist[-1])+1+safetystep
          #print "State 10 ",lo, hi 
          state = 0 #for last time
      elif state < 0:
        print state, newlist
    return newlist, statlist
  #
  #body of convergentBreaks_Inner
  #        
  if trace:
      if TraceFile !="":
        tf=open(TraceFile,"w")
      else:
        tf = None

  rt0=recursetest.recurse(testdata, controldata, datayears, model, smooth=False, trim = 0, pr=screenpr, anom=False, withshifts=True)
  
  oyears = np.sort(rt0.breakyears().keys()) #initial set of breakpoints 
  byears = [datayears[0]]
  byears.extend(oyears.tolist())
  #print "BYEARS", byears
  #byears=np.insert(byears, 0, datayears[0])
  
  if trace:
    print >>tf,model, "trace=",trace, "shallow=",shallow
    for i in range(len(testdata)):
      print >>tf,i, testdata[i], controldata[i], datayears[i], aicControl[i]
  byears.append(datayears[-1]) 
  if trace: print >>tf,"BYEARS", byears

  initialBreaks= copy.copy(byears)
  #low=0
#  crit30 = bivariate.critTi(0.01, 30) #JHR 16sep14 this gets over written further down
  #print "initially", np.sort(oyears.keys()), crit30
  solutions={}
  aics={}
  done = len(byears) <= 2
 
  
  statlist = []
  newbreaks = copy.copy(byears) 
  
  while not done:
    newbreaks, statlist=merge_pass(byears, testdata, controldata, datayears)
    #save the solutions so we can select the best in case we go too far
    key=str(newbreaks)
    if key in solutions:
      done = True
    else:
      try:
        aics[key]=AIC.AIC(testdata, datayears, datayears, newbreaks,breaksonly=True ).value()
      except Exception as e:
        if trace: print "AIC key ",key, "cannot be evaluated", str(e)
        aics[key] = 10.0e10
        pass
      solutions[key]=(initialBreaks, newbreaks, statlist, aics[key])
    if len(solutions) >2:
      done = True

    if shallow: #halt convergence is not required
      done = True
    else:
      byears=newbreaks
  #now return either the best solution or the convergent solution
  if not terminalsolution:
    aic=None
    for key in solutions.keys():
      if trace: print >>tf,key, AIC.AIC(testdata, datayears, datayears, initialBreaks,breaksonly=True ).value()
      if aic==None or aics[key][0] < aic[0]:
        (initialBreaks, newbreaks, statlist, aic) = solutions[key]
        print "prefer ", key, aic, "Initial",AIC.AIC(testdata, datayears, datayears, initialBreaks,breaksonly=True ).value()
        if trace: print >>tf,"prefer ", key, aic, "Initial",AIC.AIC(testdata, datayears, datayears, initialBreaks,breaksonly=True ).value()
  print "\nReturning",initialBreaks, "->", newbreaks, statlist
  if trace: print >>tf,"Returning",initialBreaks, "->", newbreaks, statlist
  return initialBreaks, newbreaks, statlist

  
controlmodes=set(["control","years","flat","trend"])
def convergentBreaks(testdata, controldata, datayears, model, mode="years", guide="Stability", pr=0.01, screenpr=0.05, trace=False, shallow=False):
  if not mode in controlmodes:
    raise CBException("keyword 'mode' should be in "+str(controlmodes))
  if guide in ["AIC", "Stability"]:
    if mode=="flat":
      ctrl=np.ones(np.shape(controldata))
      return convergentBreaks_Inner(testdata, controldata, datayears, ctrl, model, pr=pr, screenpr=screenpr, trace=trace, shallow=shallow)
    elif mode=="years":
      return convergentBreaks_Inner(testdata, controldata, datayears, datayears, model, pr=pr, screenpr=screenpr, trace=trace, shallow=shallow)
    elif mode=="control":
      return convergentBreaks_Inner(testdata, controldata, datayears, controldata, model, pr=pr, screenpr=screenpr, trace=trace, shallow=shallow)
    else:
      rstats = regress.analysed_regress(testdata, controldata)
      yh, r = regress.residuals(testdata, controldata, rstats)
      return convergentBreaks_Inner(testdata, controldata, datayears, yh, mode, pr=pr, screenpr=screenpr, trace=trace, shallow=shallow)

if __name__ == "__main__":
  import datetime
  import os
  
  timestring="%d-%d-%d-%d-%d" % datetime.datetime.timetuple(datetime.datetime.utcnow())[0:5]
  prewhiten=False
  fn=os.environ["HOMEPATH"]+"\\Documents\\ReferenceData\\Rogers_Analysis_17Jul2014\\Data 4 Jim\\GISSTEMPto6-2013\\GISSTEMPto6-2013b.csv"
  fn=os.environ["HOMEPATH"]+"\\Documents\\abrupt\\SpatialBPs\\Qld_191001-201312.csv"
  fn=os.environ["HOMEPATH"]+"\\Documents\\ReferenceData\\Rogers_Analysis_17Jul2014\\Data 4 Jim\\GISSTEMPto6-2013\\GISSTEMPto6-2013b.csv"
  fn="C:\Users\s4493222\Documents\CourseWork(ROP8001&2)\Confirmation\IllustrativeData.csv"
  fn=os.environ["HOMEPATH"]+"\\Documents\\ReferenceData\\GISS\\gisstemp\\tabledata_v3\\GLB.Ts+dSST.csv"
  TraceFile=fn+".trace"
  trace=True
  #fn=os.environ["HOMEPATH"]+"\\Documents\\ReferenceData\\BOM\\rutherglen_82039.csv"
  #fn = tests17Aug.fn
  data=np.genfromtxt(fn,delimiter=",",names=True,filling_values =np.NaN, skiprows=0) #fill with NaNs so must use their bounds
  ys=data["JD"]
  #ys=data["DeWobble"]
#  ys=data["24SEQU"]
  #ys=data["Mint"]
  Years=data["Year"]
  #xs=data["Maxt"]
  xs = np.array([random.random() for i in range(len(ys))])
  if prewhiten:
    import STARS
    import whitening
    ys=whitening.prewhiten(ys,STARS.AlphaEst(ys, 15, option="optIPN4", returnmsgs=False))
  print xs, ys
  print convergentBreaks(ys, xs, Years, "64N90N", mode="control", guide="AIC",trace=trace)
  
  bv=bivariate.bivariate(ys, xs, anomalise=False)
  print bv.maxTi(), bv.maxIndexTi()
  print bv.allPoints(0.05)
