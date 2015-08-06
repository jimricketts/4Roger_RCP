# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 17:26:51 2015

@author: s4493222
"""

SVNRevision="$Revision"

#[jhr579@raijin4 GCM]$ for d in *;do dd=`ls $d`;for ddd in $dd;do echo "$d:$ddd";done

modeltable={"bcc-csm1-1":"BCC",
"bcc-csm1-1-m":"BCC",
"BNU-ESM":"BNU",
"KCM1-2-2":"CAU-GEOMAR",
"CanAM4":"CCCMA",
"CanCM4":"CCCMA",
"CanESM2":"CCCMA",
"CMCC-CESM":"CMCC",
"CMCC-CM":"CMCC",
"CMCC-CMS":"CMCC",
"CNRM-CM5":"CNRM",
"CNRM-CM5-2":"CNRM-CERFACS",
"CFSv2-2011":"COLA-CFS",
#"CSIRO-Mk3-6-0":"CSIRO",
"ACCESS1-0":"CSIRO-BOM",
"ACCESS1-3":"CSIRO-BOM",
"CSIRO-Mk3-6-0":"CSIRO-QCCCE",
"FIO-ESM":"FIO",
"COSMOS-ASO":"FUB",
"EC-EARTH":"ICHEC",
"EC-EARTH-2-2":"ICHEC",
"inmcm4":"INM",
"IPSL-CM5A-LR":"IPSL",
"IPSL-CM5A-MR":"IPSL",
"IPSL-CM5B-LR":"IPSL",
"FGOALS-g2":"LASG-CESS",
"FGOALS_g2":"LASG-CESS",
"FGOALS-gl":"LASG-IAP",
"FGOALS-s2":"LASG-IAP",
"MIROC-ESM":"MIROC",
"MIROC-ESM-CHEM":"MIROC",
"MIROC4h":"MIROC",
"MIROC5":"MIROC",
"HadCM3":"MOHC",
"HadGEM2-A":"MOHC",
"HadGEM2-CC":"MOHC",
"HadGEM2-ES":"MOHC",
"MPI-ESM-LR":"MPI-M",
"MPI-ESM-MR":"MPI-M",
"MPI-ESM-P":"MPI-M",
"MRI-AGCM3-2H":"MRI",
"MRI-AGCM3-2S":"MRI",
"MRI-CGCM3":"MRI",
"MRI-ESM1":"MRI",
"GISS-E2-H":"NASA-GISS",
"GISS-E2-H-CC":"NASA-GISS",
"GISS-E2-R":"NASA-GISS",
"GISS-E2-R-CC":"NASA-GISS",
"GEOS-5":"NASA-GMAO",
"CCSM4":"NCAR",
"NorESM1-M":"NCC",
"NorESM1-ME":"NCC",
"HadGEM2-AO":"NIMR-KMA",
"GFDL-CM2p1":"NOAA-GFDL",
"GFDL-CM3":"NOAA-GFDL",
"GFDL-ESM2G":"NOAA-GFDL",
"GFDL-ESM2M":"NOAA-GFDL",
"GFDL-HIRAM-C180":"NOAA-GFDL",
"GFDL-HIRAM-C360":"NOAA-GFDL",
"CFSv2-2011":"NOAA-NCEP",
"CESM1-BGC":"NSF-DOE-NCAR",
"CESM1-CAM5":"NSF-DOE-NCAR",
"CESM1-CAM5-1-FV2":"NSF-DOE-NCAR",
"CESM1-FASTCHEM":"NSF-DOE-NCAR",
"CESM1-WACCM":"NSF-DOE-NCAR",
"CSIRO-Mk3L-1-2":"UNSW"}

groupmodels={}
for k in modeltable:
  if not modeltable[k] in groupmodels:
    groupmodels[modeltable[k]] = []
  groupmodels[modeltable[k]].append(k)
'''
groupmodels={'ICHEC': ['EC-EARTH', 'EC-EARTH-2-2'], 
'LASG-CESS': ['FGOALS-g2', 'FGOALS_g2'], 
'BCC': ['bcc-csm1-1-m', 'bcc-csm1-1'], 
'NSF-DOE-NCAR': ['CESM1-CAM5-1-FV2', 
'CESM1-CAM5', 'CESM1-BGC', 'CESM1-FASTCHEM', 'CESM1-WACCM'], 'FUB': ['COSMOS-ASO'], 'MOHC': ['HadCM3', 'HadGEM2-CC', 'HadGEM2-A', 'HadGEM2-ES'], 'IPSL': ['IPSL-CM5B-LR', 'IPSL-CM5A-MR', 'IPSL-CM5A-LR'], 'CNRM': ['CNRM-CM5'], 'CMCC': ['CMCC-CMS', 'CMCC-CM', 'CMCC-CESM'], 'CAU-GEOMAR': ['KCM1-2-2'], 'CSIRO-BOM': ['ACCESS1-0', 'ACCESS1-3'], 'MPI-M': ['MPI-ESM-P', 'MPI-ESM-MR', 'MPI-ESM-LR'], 'NCAR': ['CCSM4'], 'NIMR-KMA': ['HadGEM2-AO'], 'CSIRO-QCCCE': ['CSIRO-Mk3-6-0'], 'CCCMA': ['CanESM2', 'CanCM4', 'CanAM4'], 'BNU': ['BNU-ESM'], 'NOAA-NCEP': ['CFSv2-2011'], 'CNRM-CERFACS': ['CNRM-CM5-2'], 'NASA-GMAO': ['GEOS-5'], 'NASA-GISS': ['GISS-E2-R-CC', 'GISS-E2-H', 'GISS-E2-R', 'GISS-E2-H-CC'], 'MIROC': ['MIROC4h', 'MIROC-ESM', 'MIROC5', 'MIROC-ESM-CHEM'], 'FIO': ['FIO-ESM'], 'NOAA-GFDL': ['GFDL-ESM2M', 'GFDL-CM3', 'GFDL-CM2p1', 'GFDL-ESM2G', 'GFDL-HIRAM-C180', 'GFDL-HIRAM-C360'], 'LASG-IAP': ['FGOALS-s2', 'FGOALS-gl'], 'INM': ['inmcm4'], 'UNSW': ['CSIRO-Mk3L-1-2'], 'NCC': ['NorESM1-M', 'NorESM1-ME'], 'MRI': ['MRI-ESM1', 'MRI-AGCM3-2H', 'MRI-CGCM3', 'MRI-AGCM3-2S']}
'''
groups=set([modeltable[k] for k in modeltable])

'''
groups=set(['ICHEC', 'LASG-CESS', 'BCC', 'NSF-DOE-NCAR', 'FUB', 'MOHC', 'IPSL', 'CNRM', 'CMCC', 'CAU-GEOMAR', 'CSIRO-BOM', 'MPI-M', 'NCAR', 'NIMR-KMA', 'CSIRO-QCCCE', 'CCCMA', 'BNU', 'NOAA-NCEP', 'CNRM-CERFACS', 'NASA-GMAO', 'NASA-GISS', 'MIROC', 'FIO', 'NOAA-GFDL', 'LASG-IAP', 'INM', 'UNSW', 'NCC', 'MRI'])
'''

centres={'ICHEC':["ie1","ireland",'EC-Earth (European Earth System Model)'],
 'LASG-CESS':["cn1","china","IAP (Institute of Atmospheric Physics, Chinese Academy of Sciences, Beijing, China) and THU (Tsinghua University)"  ],
 'BCC':["cn2","china","Beijing Climate Center(BCC),China Meteorological Administration,China"],  
 'NSF-DOE-NCAR':["us1","USA","NSF/DOE NCAR (National Center for Atmospheric Research) Boulder, CO, USA" ],                 
 'FUB':["de1","germany","FUB IfM (Freie Universitaet Berlin, Institute for Meteorology, Berlin, Germany)"], 
 'MOHC':["uk1","united kingdom","a.schurer@ed.ac.uk, The King\'s Buildings, West Mains Road, Edinburgh" ], 
 'IPSL':["fr1","france","Model documentation and further reference available here : http://icmc.ipsl.fr"],
 'CNRM':["fr2","france","CNRM (Centre National de Recherches Meteorologiques, Meteo-France, Toulouse,France) and CERFACS (Centre Europeen de Recherches et de Formation Avancee en Calcul Scientifique, Toulouse, France)" ],
 'CMCC':["it1" ,"italy" , "CMCC - Centro Euro-Mediterraneo per i Cambiamenti" ], 
 'CAU-GEOMAR':["de2", "germany", "CAU (Institut fuer Geowissenschaften), Kiel, Germany" ], 
 'CSIRO-BOM':["au1","australia","CSIRO (Commonwealth Scientific and Industrial Research Organisation, Australia), and BOM (Bureau of Meteorology, Australia)" ], 
 'MPI-M':["de3","germany","Max Planck Institute for Meteorology" ], 
 'NCAR':["us2","USA","NCAR (National Center for Atmospheric Research) Boulder, CO, USA" ] ,
 'NIMR-KMA':["kr1","south korea","NIMR (National Institute of Meteorological Research, Seoul, South Korea)"], 
 'CSIRO-QCCCE':["au2","australia","Australian Commonwealth Scientific and Industrial Research Organization (CSIRO) Marine and Atmospheric Research (Melbourne, Australia) in collaboration with the Queensland Climate Change Centre of Excellence (QCCCE) (Brisbane, Australia)"], 
 'CCCMA':["ca1","canada","CCCma (Canadian Centre for Climate Modelling and Analysis, Victoria, BC, Canada)" ], 
 'BNU':["cn3","china""GCESS,BNU,Beijing,China" ], 
 'NOAA-NCEP':["us3","USA","NCEP (National Centers for Environmental Prediction, Camp Springs, MD)" ], 
 'CNRM-CERFACS':["fr3","france","CNRM (Centre National de Recherches Meteorologiques, Meteo-France, Toulouse, France) and CERFACS (Centre Europeen de Recherches et de Formation Avancee en Calcul Scientifique, Toulouse, France)"], 
 'NASA-GMAO':["us4","USA","Global Modeling and Assimilation Office, NASA Goddard Space Flight Center, Greenbelt, MD 20771"], 
 'NASA-GISS':["us5","USA","NASA/GISS (Goddard Institute for Space Studies) New York, NY"], 
 'MIROC':["jp1","japan","AORI (Atmosphere and Ocean Research Institute, The University of Tokyo, Chiba, Japan), NIES (National Institute for Environmental Studies, Ibaraki, Japan), JAMSTEC (Japan Agency for Marine-Earth Science and Technology, Kanagawa, Japan)"], 
 'FIO':["cn4","china","FIO(The First Institution of Oceanography,SOA,Qingdao,China)" ], 
 'NOAA-GFDL':["us6","USA","NOAA GFDL(201 Forrestal Rd, Princeton, NJ, 08540)" ], 
 'LASG-IAP':["cn5","china","IAP(Institute of Atmospheric Physics),CAS(Chinese Academy of Sciences),Beijing,China"], 
 'INM':["ru1","russia","INM (Institute for Numerical Mathematics,  Moscow, Russia)"], 
 'UNSW':["au3","australia","UNSW (University of New South Wales, Sydney, Australia)"], 
 'NCC':["no1","norway","Norwegian Climate Centre"] , 
 'MRI':["jp2","japan","MRI (Meteorological Research Institute, Tsukuba, Japan)"]
 }

countries=['au',  'ca', 'cn', 'de','fr', 'ie', 'it', 'jp', 'kr','no', 'ru', 'uk','us']

'''
centre_codes={}
for k in centres:
  centre_codes[centres[k][0]]=k
'''
centre_codes={'de3': 'MPI-M', 
'fr1': 'IPSL', 
'ru1': 'INM', 
'no1': 'NCC', 'au1': 'CSIRO-BOM', 'us6': 'NOAA-GFDL', 'de1': 'FUB', 'de2': 'CAU-GEOMAR', 'kr1': 'NIMR-KMA', 'au3': 'UNSW', 'us4': 'NASA-GMAO', 'jp1': 'MIROC', 'jp2': 'MRI', 'it1': 'CMCC', 'cn4': 'FIO', 'cn5': 'LASG-IAP', 'cn2': 'BCC', 'cn3': 'BNU', 'cn1': 'LASG-CESS', 'ie1': 'ICHEC', 'uk1': 'MOHC', 'fr3': 'CNRM-CERFACS', 'fr2': 'CNRM', 'us5': 'NASA-GISS', 'au2': 'CSIRO-QCCCE', 'us3': 'NOAA-NCEP', 'us2': 'NCAR', 'us1': 'NSF-DOE-NCAR', 'ca1': 'CCCMA'}
 
'''
print groups
print modeltable.keys()
a=[centres[k][0] for k in centres]
b=set(a)
print len(a), len(b), a, centre_codes,groupmodels
'''
import glob
import os

def codes(model, replicate):
  #for a given mode,find its group, e.g. given EC-EARTH, return ICHEC
  grp = modeltable[model]
  #for a group, find the group alpahnum code
  alphanumcode = centres[grp][0]
  #for that code, assign a country number 
  cnumber = 1 + int(countries.index(alphanumcode[:2]))
  #and groupincointryno
  groupincountryno=int(alphanumcode[2:]) 
  #now the modelingroupnumber
  modelingroupnumber=groupmodels[grp].index(model)+1
  #groupno is cnumber*10+groupincountryno
  groupno=cnumber*10+groupincountryno
  #now the replicate
  r=replicate[1:].split('i')
  i=r[1].split('p')
  p=i[1]
  i=i[0]
  r=r[0]
  rip="%04d" % (100*int(r)+10*int(i)+int(p),)
    
  return grp,alphanumcode,"%03d"%(groupno,),model, "%02d"% (modelingroupnumber,), "%03d%02d"%(groupno,modelingroupnumber), replicate, rip
  

def matchfn(fn):
  pos=-1
  models = modeltable.keys()
  m=0
  while m <len(models) and fn.find(models[m]) <0:
    m+=1
  if m < len(models):
    modl=models[m]
    return codes(modl,os.path.basename(fn).split("_")[0])    
  else:
    return None

if __name__ == "__main__":
  fns=glob.glob("c:\\users\\s4493222\\documents\\abrupt\\c_test\\r*.dat")
  print fns
  
  for fn in fns:
    code=matchfn(fn)
    if code != None:
      print fn, code
