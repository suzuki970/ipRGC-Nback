#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 16:52:13 2020

@author: yuta
"""

import numpy as np
import glob
import os
import json

def au2mm(visDist):

    fileName = glob.glob(os.path.dirname(os.path.realpath(__file__)) + "/au2mm/data/Eyelink/" + str(visDist) + "/asc/*")
    fileName.sort()
    
    datHash = {}
    for n in [6,8]:
        datHash[str(n)] = []
        
    for file in fileName:
        
        f = open(os.path.join(str(file)))
    
        dat=[]
        for line in f.readlines():
            dat.append(line.split())
            
        f.close()
        
        tmp = []
        for i,line in enumerate(dat):
            if i > 30:
                if line[0].isdecimal():
                    tmp.append(float(line[3]))
       
        if os.name == "nt":
            datHash[file.split('\\')[-1][-7:-6]] .append(tmp[5000:10000])
        else:
            datHash[file.split('/')[-1][-7:-6]] .append(tmp[5000:10000])
                
    
    mmName = list(datHash.keys())
    
    dat = []
    for mm in mmName:
        for d in datHash[mm]:
            dat.append( int(mm) / np.sqrt(np.mean(d) ))
        
    return np.mean(dat)

def getmm(eyeTracker = "VSHD"):

    fileName = glob.glob(os.path.dirname(os.path.realpath(__file__)) + "/au2mm/data/" + eyeTracker + "/*.txt")
    # fileName = glob.glob("./au2mm/data/" + eyeTracker + "/*.txt")
    fileName.sort()
    
    datHash = {}
    for n in [3,7]:
        datHash[str(n)+"mm"] = []
    
    for iFile in np.arange(len(fileName)):
        f = open(fileName[iFile],encoding = "shift_jis")
        dat=[]
        for line in f.readlines():
            dat.append(line.split())        
        f.close()
            
        tmp = []
        for i,line in enumerate(dat):
            if (i > 44) and ('MovieFrame' not in line) and ("eyeA:gazeNudge" not in line) and ("eyeB:gazeNudge" not in line) and ('END' not in line) :
                
                tmp.append([float(line[8]),   # width
                            float(line[9])])  # height
        
        if os.name == "nt":                 
            datHash[fileName[iFile].split('\\')[-1][0]+"mm"] = np.array(tmp).mean(axis=0)
        else:
            datHash[fileName[iFile].split('/')[-1][0]+"mm"] = np.array(tmp).mean(axis=0)
        
    dat = []
    for mm in list(datHash.keys()):
        # for d in datHash[mm]:
        dat.append( int(mm[0]) / datHash[mm] )
    
    return np.array(dat).mean(axis=0)


    #%%
    
    # fileName = glob.glob(os.path.dirname(os.path.realpath(__file__)) + "/au2mm/data/" + eyeTracker + "/2mm_LeftEye.mp4_eyeRec_table.json")
    # # fileName = glob.glob("./au2mm/data/" + eyeTracker + "/2mm_LeftEye.mp4_eyeRec_table.json")
    # fileName.sort()
    # with open(fileName[0], 'r') as f:
    #     dat = json.load(f)
    
    # datHash={"3mm":[]}
    # for d in dat:
    #     datHash["3mm"].append( 3/ d["size"][0])
    
    # return np.array(datHash["3mm"]).mean(axis=0)

