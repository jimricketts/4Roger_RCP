# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 16:36:19 2015

@author: s4493222
"""
SVNRevision="$Revision: 307 $"
from HADCRUT4 import HADCRUT4

class COWTAN_WAY(HADCRUT4):
  def __init__(self, filename,annual=True):
    super(COWTAN_WAY,self).__init__(filename, annual=annual)
