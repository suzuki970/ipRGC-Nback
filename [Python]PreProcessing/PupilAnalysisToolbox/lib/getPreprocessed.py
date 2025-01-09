#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 09:42:22 2022

@author: yutasuzuki
"""

import numpy as np
import matplotlib.pyplot as plt
import json
from asc2array_cls import asc2array_cls
import glob
import os
import datetime
import pandas as pd
from pre_processing_cls import pre_processing,rejectDat,rejectedByOutlier
import matplotlib
import seaborn as sns        
from LightSource import LightSource,showSpectra,getColor
 
from scipy import interpolate
from pre_processing_cls import getNearestValue

import math
from pixel_size import pixel2angle
from multiprocessing import Pool

dt_now = datetime.datetime.now()
date = dt_now.strftime("%Y%m%d")

rootFolder = "../../experimental_scripts/3A_GlareIllusion/results/"
   
folderList=[]

for filename in os.listdir(rootFolder):
    if os.path.isdir(os.path.join(rootFolder, filename)): 
        folderList.append(filename)
folderList.sort()

cfg={   
     "TIME_START":-10,
     "TIME_END":4,
     "WID_ANALYSIS":4,
     "WID_BP_ANALYSIS":2,
     "useEye":1,
     "WID_FILTER":[],
     "RESAMPLING_RATE":1000,
     "mmFlag":False,
     "normFlag":False,
     "s_trg":[],
     "visualization":False,
     "MS":False,
     "DOT_PITCH":0.33,   
     "VISUAL_DISTANCE":100,
     "acceptMSRange":4.5,
     "MIN_DURATION":0.1,
     "SCREEN_RES":[640*2, 512*2],
     "rejectFlag":[]
     }

center = np.array(cfg['SCREEN_RES'])/2
t2a = asc2array_cls(cfg)
  
if not cfg['mmFlag'] and not cfg['normFlag']:
    unitName = "_au"
elif cfg['mmFlag']:
    unitName = "_mm"
else:
    unitName = "_norm" 
    

def run(iSub):
    
    #%%
    datHash={"PDR":[],
             "PDR_all":[],
             # "Baseline":[],
             "gazeX":[],
             "gazeY":[],
             "ipRGC":[],
             "Run":[],
             "sub":[],
             "light":[],
             "pattern":[]
              }

    condition_list = {"light":[],
                      "pattern":[]}

    datHashAll = {"PDR":[],
                  "ipRGC":[]
                  }
    cfg["rejectFlag"] = []

    #% ------------------ load json file (eye data) -------------------
    print("Processing --> " + folderList[iSub][-3:] + "...")
    
    ascFileName = glob.glob(os.path.join(rootFolder + "/" + folderList[iSub] + "/*.asc"))    
    ascFileName.sort()
    
    for iRun,ascFile in enumerate(ascFileName):
        
        #% ------------------ load json file (eye data) -------------------
        dat = t2a.dataExtraction(ascFile)
    
        #% ------------------ load json file (eye data) -------------------
        eyeData,events,initialTimeVal,fs = t2a.dataParse(dat)
        ave,sigma = t2a.getAve(eyeData["pupilData"])
        
        #% ------------------ load json file (eye data) -------------------
        eyeData = t2a.blinkInterp(eyeData)
        
        if cfg["normFlag"]:
            # pupilData = t2a.pupilNorm(eyeData["pupilData"], normDat[normDat["sub"]==iSub]["ave"][iSub], normDat[normDat["sub"]==iSub]["sigma"][iSub]).reshape(-1)
            pupilData = t2a.pupilNorm(eyeData["pupilData"], ave, sigma).reshape(-1)
        else:
            
            pupilData = np.mean(eyeData["pupilData"],axis=0)
            pupilData = pupilData.reshape(-1)
    
        # pupilData = eyeData["pupilData"]
        gazeX = np.mean(eyeData["gazeX"],axis=0)
        gazeY = np.mean(eyeData["gazeY"],axis=0)
        
        cfg["rejectFlag"].append(eyeData["rejectFlag"])
       
        cfg["SAMPLING_RATE"] = int(fs)
        
        # MSC_table = {
        #     "timestamp":np.arange(eyeData["gazeX"].shape[1]),
        #     "horizon_left":eyeData["gazeX"][0,] - cfg["SCREEN_RES"][0]/2,
        #     "vertical_left": eyeData["gazeY"][0,] - cfg["SCREEN_RES"][1]/2,
        #     "horizon_right":eyeData["gazeX"][1,] - cfg["SCREEN_RES"][0]/2,
        #     "vertical_right":eyeData["gazeY"][1,] - cfg["SCREEN_RES"][1]/2,
        #     "blink": eyeData["gazeY"][0,].copy()
        #     }
       
        # MSC_table["blink"][MSC_table["blink"]>0] = 1
        
        # MSC_table["timestamp"] = np.arange(eyeData["gazeX"].shape[1])
        # MSC_table["horizon_left"] = eyeData["gazeX"][0,]
        # MSC_table["vertical_left"] = eyeData["gazeY"][0,]
        # MSC_table["horizon_right"] = eyeData["gazeX"][1,]
        # MSC_table["vertical_right"] = eyeData["gazeY"][1,]
        # MSC_table["blink"] = eyeData["gazeY"][0,]
        # MSC_table.to_json("MSC_table.json")
        
        # for mm in list(MSC_table.keys()):
        #     if not isinstance(MSC_table[mm],list):
        #         MSC_table[mm] = MSC_table[mm].tolist()
        #     if mm != "timestamp" and mm != "blink":
        #         MSC_table[mm] = pixel2angle(cfg["DOT_PITCH"],MSC_table[mm],cfg["VISUAL_DISTANCE"]).tolist()
        
        # with open(os.path.join("MSC_table.json"),"w") as f:
        #     json.dump(MSC_table,f)
   
        #% ------------------ reject if no answer -------------------
        seq_onset  = [[int((int(e[0])- initialTimeVal) * (fs/1000)),e[1]] for e in events["MSG"] if e[1] == "stim_onset"]       
        
        for c in list(condition_list.keys()):
            condition_list[c] = [int(e[2]) for e in events["MSG"] if "Condition_" + c in e[1]]           
       
        #% ------------------ data extraction ------------------   
        for s in seq_onset:
             datHash["sub"].append(int(folderList[iSub][1:3]))       
                    
        for s in seq_onset:
            datHash["PDR"].append(pupilData[np.arange(int(s[0]+fs*cfg["TIME_START"]),int(s[0]+fs*cfg["TIME_END"]))].tolist())
       
        #% ------------------ sub number ----------------------- 
        # datHash["sub"] = np.r_[datHash["sub"],np.ones(len(condition))*int(folderList[iSub][1:3])]      
        datHash["Run"] = np.r_[datHash["Run"],np.ones(len(seq_onset))*(iRun+1)]
        
        for c in list(condition_list.keys()):
            datHash[c] = np.r_[datHash[c],np.array(condition_list[c])]
    
    for mm in list(datHash.keys()):
        if not isinstance(datHash[mm],list):
            datHash[mm] = datHash[mm].tolist()
              
    return datHash
