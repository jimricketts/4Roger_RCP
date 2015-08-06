# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 09:43:20 2015

@author: James
"""
SVNRevision="$Revision: 308 $"
# for d in */.;do for f in $d/*.trace;do grep Final\ pass\ found $f;if [ $? == 0 ];then grep Return $f;echo $f;fi;done;done>x.x
import statbreaks

with open("x.x","r") as f:
  state = 0
  lines = f.readlines()
  for line in lines:

    if line[:5] =="Final":
      state = 1
    elif line[:len("Returning")] == "Returning":
      state = 2
    else:
      state = 3
      
    if state == 1:
      #find 1988 1996] 
      words = line.split()
      word=words[0]
      while word[0] == "[" or word[-1]!= "]":
        word=words.pop(0)
      try:
        if word[-2] == ".":
          passyr=float(word[:-2].strip())
        else:
          passyr=float(word[:-1].strip())
      except:
        passyr=word[:-2].strip()
        pass
    if state == 2:
      #find the last date should be returned
      words = line.split()
      word=words[0]
      while word[0] == "[" or word[-1]!= "]":
        word=words.pop(0)
      try:
        if word[-2] == ".":
          finalyr=float(word[:-2].strip())
        else:
          finalyr=float(word[:-1].strip())
      except:
        finalyr=word[:-2].strip()
        pass

    if state==3:
      filename=line
      state = 4
    if state == 4:
      if passyr!= finalyr:
        sf=statbreaks.brkrpt(filename[:-1])
        brks=sf.breaks()
        if passyr > brks[-2]:
          print passyr== finalyr,passyr, finalyr, brks[-2:], filename[:-1]
      state=0
      
