# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 16:43:22 2015

@author: s4493222
"""

SVNRevision="$Revision$"
import cmip5_info
cmip3names=[["bccr_bcm2_0.sresa1b.run1","BCC-A1-r1","BCCA1r1"],
["bcc_cm1.sresa2.run1","BCC-A2-r1","BCCA2r1"],
["cccma_cgcm3_1.sresa1b.E1_2_3","CCCMA-A1-E","CCCMAA1E"],
["cccma_cgcm3_1.sresa1b.run1","CCCMA-A1-r1","CCCMAA1r1"],
["cccma_cgcm3_1.sresa1b.run2","CCCMA-A1-r2","CCCMAA1r2"],
["cccma_cgcm3_1.sresa1b.run3","CCCMA-A1-r3","CCCMAA1r3"],
["cccma_cgcm3_1.sresa1b.run4","CCCMA-A1-r4","CCCMAA1r4"],
["cccma_cgcm3_1.sresa1b.run5","CCCMA-A1-r5","CCCMAA1r5"],
["cccma_cgcm3_1.sresa2.run1","CCCMA-A2-r1","CCCMAA2r1"],
["cccma_cgcm3_1_t63.sresa1b.run1","CCCMA-t63-A1-r1","CCCMAt63A1r1"],
["cnrm_cm3.sresa1b.run1","CNRM-A1-r1","CNRMA1r1"],
["cnrm_cm3.sresa2.run1","CNRM-A2-r1","CNRMA2r1"],
["csiro_mk3_0.sresa1b.run1","CSIRO3-A1-r1","CSIRO3A1r1"],
["csiro_mk3_0.sresa2.run1","CSIRO3-A2-r1","CSIRO3A2r1"],
["csiro_mk3_5.sresa1b.run1","CSIRO35-A1-r1","CSIRO35A1r1"],
["csiro_mk3_5.sresa2.run1","CSIRO35-A2-r1","CSIRO35A2r1"],
["gfdl_cm2_0.sresa1b.run1","GFDL2-A1-r1","GFDL2A1r1"],
["gfdl_cm2_0.sresa2.run1","GFDL2-A2-r1","GFDL2A2r1"],
["gfdl_cm2_1.sresa1b.run1","GFDL21-A1-r1","GFDL21A1r1"],
["gfdl_cm2_1.sresa2.run1","GFDL21-A2-r1","GFDL21A2r1"],
["giss_aom.sresa1b.E1_2","GISS-A1-E","GISSA1E"],
["giss_aom.sresa1b.run1","GISS-A1-r2","GISSA1r2"],
["giss_aom.sresa1b.run2","GISS-A1-r1","GISSA1r1"],
["giss_model_e_h.sresa1b.E1_2_3","GISSE-A1-E","GISSEA1E"],
["giss_model_e_h.sresa1b.run1","GISSE-A1-r1","GISSEA1r1"],
["giss_model_e_h.sresa1b.run2","GISSE-A1-r2","GISSEA1r2"],
["giss_model_e_h.sresa1b.run3","GISSE-A1-r3","GISSEA1r3"],
["giss_model_e_r.sresa1b.E2_4","GISSE-A1-r4","GISSEA1r4"], ##INCORRECT CODING
["giss_model_e_r.sresa1b.run1","GISSER-A1-r1","GISSERA1r1"],
["giss_model_e_r.sresa1b.run2","GISSER-A1-r2","GISSERA1r2"],
["giss_model_e_r.sresa1b.run3","GISSER-A1-r3","GISSERA1r3"],
["giss_model_e_r.sresa1b.run4","GISSER-A1-r4","GISSERA1r4"],
["giss_model_e_r.sresa2.run1","GISSER-A2-r1","GISSERA2r1"],
["iap_fgoals1_0_g.sresa1b.E1_2_3","IAP-A1-E","IAPA1E"],
["iap_fgoals1_0_g.sresa1b.run1","IAP-A1-r1","IAPA1r1"],
["iap_fgoals1_0_g.sresa1b.run2","IAP-A1-r2","IAPA1r2"],
["iap_fgoals1_0_g.sresa1b.run3","IAP-A1-r3","IAPA1r3"],
["ingv_echam4.sresa1b.run1","INGV-A1-r1","INGVA1r1"],
["ingv_echam4.sresa2.run1","INGV-A2-r1","INGVA2r1"],
["inmcm3_0.sresa1b.run1","INMCM3-A1-r1","INMCM3A1r1"],
["inmcm3_0.sresa2.run1","INMCM3-A2-r1","INMCM3A2r1"],
["ipsl_cm4.sresa1b.run1","IPSL-A1-r1","IPSLA1r1"],
["ipsl_cm4.sresa2.run1","IPSL-A2-r1","IPSLA2r1"],
["miroc3_2_hires.sresa1b.run1","MIROChi-A1-R1","MIROChiA1R1"],
["miroc3_2_medres.sresa1b.E2_3","MIROC-A1-E","MIROCA1E"],
["miroc3_2_medres.sresa1b.run1","MIROC-A1-r1","MIROCA1r1"],
["miroc3_2_medres.sresa1b.run2","MIROC-A1-r2","MIROCA1r2"],
["miroc3_2_medres.sresa1b.run3","MIROC-A1-r3","MIROCA1r3"],
["miroc3_2_medres.sresa2.E1_2_3","MIROC-A2-E","MIROCA2E"],
["miroc3_2_medres.sresa2.run1","MIROC-A2-r1","MIROCA2r1"],
["miroc3_2_medres.sresa2.run2","MIROC-A2-r2","MIROCA2r2"],
["miroc3_2_medres.sresa2.run3","MIROC-A2-r3","MIROCA2r3"],
["miub_echo_g.sresa1b.E1_2_3","MIUB-A1-E","MIUBA1E"],
["miub_echo_g.sresa1b.run1","MIUB-A1-r1","MIUBA1r1"],
["miub_echo_g.sresa1b.run2","MIUB-A1-r2","MIUBA1r2"],
["miub_echo_g.sresa1b.run3","MIUB-A1-r3","MIUBA1r3"],
["miub_echo_g.sresa2.E1_2_3","MIUB-A2-E","MIUBA2E"],
["miub_echo_g.sresa2.run1","MIUB-A2-r1","MIUBA2r1"],
["miub_echo_g.sresa2.run2","MIUB-A2-r2","MIUBA2r2"],
["miub_echo_g.sresa2.run3","MIUB-A2-r3","MIUBA2r3"],
["mpi_echam5.sresa1b.E1_2_3","MPI-A1-E","MPIA1E"],
["mpi_echam5.sresa1b.run1","MPI-A1-r1","MPIA1r1"],
["mpi_echam5.sresa1b.run2","MPI-A1-r2","MPIA1r2"],
["mpi_echam5.sresa1b.run3","MPI-A1-r3","MPIA1r3"],
["mpi_echam5.sresa1b.run4","MPI-A1-r4","MPIA1r4"],
["mpi_echam5.sresa2.E1_2_3","MPI-A2-E","MPIA2E"],
["mpi_echam5.sresa2.run1","MPI-A2-r1","MPIA2r1"],
["mpi_echam5.sresa2.run2","MPI-A2-r2","MPIA2r2"],
["mpi_echam5.sresa2.run3","MPI-A2-r3","MPIA2r3"],
["mri_cgcm2_3_2a.sresa1b.E1_2_3_4_5","MRI-A1-E","MRIA1E"],
["mri_cgcm2_3_2a.sresa1b.run1","MRI-A1-r1","MRIA1r1"],
["mri_cgcm2_3_2a.sresa1b.run2","MRI-A1-r2","MRIA1r2"],
["mri_cgcm2_3_2a.sresa1b.run3","MRI-A1-r3","MRIA1r3"],
["mri_cgcm2_3_2a.sresa1b.run4","MRI-A1-r4","MRIA1r4"],
["mri_cgcm2_3_2a.sresa1b.run5","MRI-A1-r5","MRIA1r5"],
["mri_cgcm2_3_2a.sresa2.E1_2_3_4_5","MRI-A2-E","MRIA2E"],
["mri_cgcm2_3_2a.sresa2.run1","MRI-A2-r1","MRIA2r1"],
["mri_cgcm2_3_2a.sresa2.run2","MRI-A2-r2","MRIA2r2"],
["mri_cgcm2_3_2a.sresa2.run3","MRI-A2-r3","MRIA2r3"],
["mri_cgcm2_3_2a.sresa2.run4","MRI-A2-r4","MRIA2r4"],
["mri_cgcm2_3_2a.sresa2.run5","MRI-A2-r5","MRIA2r5"],
["ncar_ccsm3_0.sresa1b.E1_2","NCAR-A1-E","NCARA1E"],
["ncar_ccsm3_0.sresa1b.run1","NCAR-A1-r1","NCARA1r1"],
["ncar_ccsm3_0.sresa1b.run2","NCAR-A1-r2","NCARA1r2"],
["ncar_ccsm3_0.sresa1b.run3","NCAR-A1-r3","NCARA1r3"],
["ncar_ccsm3_0.sresa1b.run4","NCAR-A1-r4","NCARA1r4"],
["ncar_ccsm3_0.sresa1b.run9","NCAR-A1-r9","NCARA1r9"],
["ncar_ccsm3_0.sresa2.run1","NCAR-A2-r1","NCARA2r1"],
["ncar_ccsm3_0.sresa2.run2","NCAR-A2-r2","NCARA2r2"],
["ncar_ccsm3_0.sresa2.run3","NCAR-A2-r3","NCARA2r3"],
["ncar_ccsm3_0.sresa2.run4","NCAR-A2-r4","NCARA2r4"],
["ncar_pcm1.sresa1b.run1","PCM-A1-r1","PCMA1r1"],
["ncar_pcm1.sresa1b.run4","PCM-A1-r4","PCMA1r4"],
["ncar_pcm1.sresa2.E2_3_4","PCM-A2-E","PCMA2E"],
["ncar_pcm1.sresa2.run1","PCM-A2-r1","PCMA2r1"],
["ncar_pcm1.sresa2.run2","PCM-A2-r2","PCMA2r2"],
["ncar_pcm1.sresa2.run3","PCM-A2-r2","PCMA2r3"], ###NB DUPLICATION IN XLSX FILE SO NUMBERS PAST HERE ARE ODD
["ncar_pcm1.sresa2.run4","PCM-A2-r2","PCMA2r4"], ###NB DUPLICATION
["ukmo_hadcm3.sresa1b.run1","HADCM3-A1-r1","HADCM3A1r1"],
["ukmo_hadcm3.sresa2.run1","HADCM3-A2-r1", "HADCM3A2r2"],
["ukmo_hadgem1.sresa1b.run1","HADGEM-A1-r1","HADGEMA1r3"],
["ukmo_hadgem1.sresa2.run1","HADGEM-A2-r1","HADGEMA2r4"]]

knownby={}
for line in cmip3names:
  knownby[line[2]] = line[:2]
  


nameset=set(['ipsl_cm4', 'miub_echo_g', 'ukmo_hadgem1', 'ncar_pcm1', 'cccma_cgcm3_1', 'gfdl_cm2_1', 'gfdl_cm2_0', 'cnrm_cm3', 'bccr_bcm2_0', 'cccma_cgcm3_1_t63', 'iap_fgoals1_0_g', 'bcc_cm1', 'inmcm3_0', 'giss_aom', 'ukmo_hadcm3', 'ncar_ccsm3_0', 'csiro_mk3_5', 'csiro_mk3_0', 'miroc3_2_medres', 'mri_cgcm2_3_2a', 'giss_model_e_h', 'miroc3_2_hires', 'ingv_echam4', 'giss_model_e_r', 'mpi_echam5'])
centres={
'BCC':["cn2","china","Beijing Climate Centre"],
'BCCR':["no2","norway","BCCR (Bjerknes Centre for Climate Research) University of Bergen, Norway (www.bjerknes.uib.no) NERSC (Nansen Environmental and Remote Sensing Center, Norway (www.nersc.no)"],
'CCCMA':["ca1","canada","CCCma (Canadian Centre for Climate Modelling and Analysis, Victoria, BC, Canada)"],
'CNRM':["fr2","france","CNRM (Centre National de Recherches Meteorologiques, Meteo-France, Toulouse, France)"],
'CSIRO':["au4","australia","CSIRO (CSIRO Atmospheric Research, Melbourne, Australia)"],
'GFDL':["us7","USA","NOAA GFDL (US Dept of Commerce / NOAA / Geophysical Fluid Dynamics Laboratory, Princeton, NJ, USA)"],         
'GISS':["us8","USA","GISS (NASA/Goddard Institute for Space Studies, New York, USA)"],
'IAP':["cn6","china","IAP(LASG, Institute of Atmospheric Physics, P.O. Box 9804, Beijing 100029, P.R. China)"],
'INGV':["it2","italy","INGV (National Institute of Geophysics and Volcanology, Bologna, Italy)"],
'INM':["ru1","russia","INM (Iistitute for Numerical Mathematics, Moscow, Russia)"],
'IPSL':["fr1","france","IPSL (Institut Pierre Simon Laplace, Paris, France)"],
'MIROC':["jp1","japan","CCSR/NIES/FRCGC (Center for Climate System Research, Tokyo, Japan / National Institute for Environmental Studies, Ibaraki, Japan / Frontier Research Center for Global Change, Kanagawa, Japan)"],
'MIUB':["de4","germany","MIUB (University of Bonn, Bonn, Germany)"],
'MPI':["de5","germany","MPI (Max Planck Institute for Meteorology,Hamburg, Germany)"],
'MRI':["jp2","japan","MRI (Meteorological Research Institute, Tsukuba, Ibaraki, Japan)"],
'NCAR':["us2","USA","NCAR (National Center for Atmospheric Research, Boulder, CO, USA)"],
'UKMO':['uk2','united kingdom',"Met Office (Exeter, Devon, EX1 3PB, UK)"]
}
'''
names=[]
for line in cmip3names:
  print line[0].split(".")
  names.append(line[0].split(".")[0])
print set(names)  

keys=[]
for line in cmip3names:
  print line[2]
  if line[2] in keys: print line[2], "DUP"
  keys.append(line[2])
print set(keys)


for line in cmip3names:
  name=line[0].upper()
  found = False
  for c in centres.keys():
    if name.find(c) >=0:
      found = True
      res=(c, name, line)
  if not found: print line[0]
print "\n--------------------------\n"   
import cmip5_info
for c in centres:
  if c in cmip5_info.centres:
    print c, "in both", centres[c][0],cmip5_info.centres[c][0]

  else:
    cocode=centres[c][0][:2]
    for c5 in cmip5_info.centres:
      if cocode == cmip5_info.centres[c5][0][:2]:
        print c,c5, cmip5_info.centres[c5]
'''        

'''
for k in centres:
  centre_codes[centres[k][0]]=k
print centre_codes
'''
centre_codes={'cn6': 'IAP', 'cn2': 'BCC', 'au4': 'CSIRO', 'ca1': 'CCCMA', 'de4': 'MIUB', 'uk2': 'UKMO', 'de5': 'MPI', 'us7': 'GFDL', 'fr2': 'CNRM', 'fr1': 'IPSL', 'us2': 'NCAR', 'j2p': 'MRI', 'jp1': 'MIROC', 'ru1': 'INM', 'it2': 'INGV', 'no2': 'BCCR', 'us8': 'GISS'}
'''
countries=set([c[:2] for c in centre_codes])
print sorted(list(countries))
'''
countries=['au', 'ca', 'cn', 'de', 'fr', 'it', 'jp', 'no', 'ru', 'uk', 'us']
'''
modeltable={}
for m in cmip3names:
  model=m[0].split('.')[0]
  #print model
  group=model.split("_")[0].upper()
  if not model in modeltable:
    modeltable[model]=group

print modeltable  

'''
modeltable={
'ipsl_cm4': 'IPSL', 
'miub_echo_g': 'MIUB', 
'ukmo_hadgem1': 'UKMO', 
'ncar_pcm1': 'NCAR', 
'cccma_cgcm3_1': 'CCCMA', 
'gfdl_cm2_1': 'GFDL', 
'gfdl_cm2_0': 'GFDL', 
'cnrm_cm3': 'CNRM', 
'bccr_bcm2_0': 'BCCR', 
'cccma_cgcm3_1_t63': 'CCCMA', 
'iap_fgoals1_0_g': 'IAP', 
'bcc_cm1': 'BCC', 
'inmcm3_0': 'INM', 
'giss_aom': 'GISS', 
'ukmo_hadcm3': 'UKMO', 
'ncar_ccsm3_0': 'NCAR', 
'csiro_mk3_5': 'CSIRO', 
'csiro_mk3_0': 'CSIRO', 
'miroc3_2_medres': 'MIROC', 
'mri_cgcm2_3_2a': 'MRI', 
'giss_model_e_h': 'GISS', 
'miroc3_2_hires': 'MIROC', 
'ingv_echam4': 'INGV', 
'giss_model_e_r': 'GISS', 
'mpi_echam5': 'MPI'
}

'''
groupmodels={}
for k in modeltable:
  if not modeltable[k] in groupmodels:
    groupmodels[modeltable[k]] = []
  groupmodels[modeltable[k]].append(k)

print groupmodels
'''

groupmodels={
'GISS': ['giss_aom', 'giss_model_e_h', 'giss_model_e_r'], 
'CCCMA': ['cccma_cgcm3_1', 'cccma_cgcm3_1_t63'], 
'MPI': ['mpi_echam5'], 
'CSIRO': ['csiro_mk3_5', 'csiro_mk3_0'], 
'BCC': ['bcc_cm1'], 
'NCAR': ['ncar_pcm1', 'ncar_ccsm3_0'], 
'BCCR': ['bccr_bcm2_0'], 
'MIROC': ['miroc3_2_medres', 'miroc3_2_hires'], 
'MRI': ['mri_cgcm2_3_2a'], 
'INM': ['inmcm3_0'], 
'UKMO': ['ukmo_hadgem1', 'ukmo_hadcm3'], 
'INGV': ['ingv_echam4'], 
'MIUB': ['miub_echo_g'], 
'CNRM': ['cnrm_cm3'], 
'IAP': ['iap_fgoals1_0_g'], 
'IPSL': ['ipsl_cm4'], 
'GFDL': ['gfdl_cm2_1', 'gfdl_cm2_0']
}
import os

def codes(model):
  modl, sres, replicate = model.split(".")
  #for a given mode,find its group, e.g. given EC-EARTH, return ICHEC
  grp = modeltable[modl]
  #for a group, find the group alpahnum code
  alphanumcode = centres[grp][0]
  #for that code, assign a country number 
  cnumber = 1 + int(cmip5_info.countries.index(alphanumcode[:2]))
  #and groupincointryno
  groupincountryno=int(alphanumcode[2:]) 
  #now the modelingroupnumber
  try:
    modelingroupnumber=groupmodels[grp].index(modl)+1
  except:
    print grp,model
    raise
  #groupno is cnumber*10+groupincountryno
  groupno=cnumber*10+groupincountryno
  #now the replicate
  if replicate[0] in ["e","E"] :
    rip = "9999"
  else: rip="%04d" % (int(replicate[3:]),)
    
  return grp,alphanumcode,"%03d"%(groupno,),modl, "%02d"% (modelingroupnumber,), "%03d%02d"%(groupno,modelingroupnumber), replicate, rip, sres[4:]
  

def matchfn(fn):
  modelcode=fn.split("_")[-1]
  initfn=knownby[modelcode]
  return codes(initfn[0])

if __name__ == "__main__":
  import glob
  fns=glob.glob("c:\\users\\s4493222\\documents\\abrupt\\c_test\\Global temps.csv*[01234567890E]")
  for fn in fns:
    print fn, matchfn(fn)
  
#c:\users\s4493222\documents\abrupt\c_test\r9i1p1_CSIRO-Mk3-6-0_RCP2.6_icmip5_tas_Amon_ens_rcp26_0-360E_-90-90N_n_026.dat 
#('CSIRO-QCCCE', 'au2',         '012',     'CSIRO-Mk3-6-0', '01',          '01201',   'r9i1p1', '1109')
#c:\users\s4493222\documents\abrupt\c_test\Global temps.csv_PCMA2r4 
#('GRP'          'country/grp'  'cygpcode' 'model'           'modnoIngrp'  'grpdlcode' 'replicate', 'pir')
#('NCAR',        'us2',         '132',      'ncar_pcm1',     '01',         '13201',    'run4',      '0004')