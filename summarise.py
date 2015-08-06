# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 15:50:26 2015

@author: s4493222
"""

import numpy as np
lines = []
#with open('GISSTEMP3.txt','r') as f:
with open('NCDC.txt','r') as f:
  lines = f.readlines()
  f.close()

mods={}
for line in lines:
  
  pos = line.find("[")
  model=line[:pos]
  if not model in mods:
    mods[model] = {}
  ylist=line[pos:-1]
  if not ylist in mods[model]:
    mods[model][ylist]=0
  mods[model][ylist] += 1

for model in mods:
  for ylist in mods[model]:
    print model, ylist, mods[model][ylist]
