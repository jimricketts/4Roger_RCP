# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 12:53:53 2015

@author: s4493222
"""

SVNRevision="$Revision: 373 $"




def doit(fn,outfn):
  with open(fn) as f:
    lines=f.readlines()
  l=0
  table={}
  for line in lines:
    l += 1
    if line[:4] != '"Seq':
      words = line.split('"')
      fields=words[2].split(",")
      datasetno=int(fields[1])
      setnoindata=int(fields[3])
      pctindata=int(fields[7])
      if not datasetno in table:
        table[datasetno] = {}
      table[datasetno][pctindata]=setnoindata
  #print table
  with open(outfn,"w") as outf:
    for line in lines:
      l += 1
      if line[:4] != '"Seq':
        words = line.split('"')
        fields=words[2].split(",")
        datasetno=int(fields[1])
        setnoindata=int(fields[3])
        pctindata=int(fields[7])
        if len(table[datasetno]) == 1:
          line2 = line[:-1]+',"Y"'
        else:
          #what is the maximimum
          maxsetnoindata=table[datasetno][max(table[datasetno].keys())]
          
          if setnoindata==maxsetnoindata:
            line2 = line[:-1]+',"Y"'
          else:
            line2 = line[:-1]+',"N"'
      else:
        line2=line[:-1]+',"Modal"'
      print >>outf,line2
    
  
if __name__ == "__main__":
  doit("piControl.txt","piControl.txt.csv")