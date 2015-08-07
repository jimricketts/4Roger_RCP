# -*- coding: utf-8 -*-
"""
Created on Thu Aug 06 16:41:26 2015

@author: s4493222
"""

SVNRevision="$Revision$"

#script to perform statistics for Roger (see ticket #24).

#1 we need to open all named csv files and create a simple data base from them.

import csv

fn="C:\\Users\\s4493222\\Documents\\abrupt\\4Roger_RCP\\RCP85_05aug.x.txt.csv"
with open(fn, 'rb') as csvfile:
  lines = csv.reader(csvfile, delimiter=',', quotechar='"')
  for line in lines:
    print line